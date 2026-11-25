"""Cell-Hub end-to-end UI test: drives the running app in a headless Chromium, module by module, on the Zenodo demo data.
Every step (tab change, input, button click, download) is recorded with its status, duration, error text and a screenshot on failure.
Usage: python ui_test.py --url http://cellhub_app:3838 --demo /demo --out /out [--modules single,cellchat,trajectory,multiple,spatial]
"""
import argparse, json, os, re, time, traceback
from playwright.sync_api import sync_playwright, TimeoutError as PWTimeout

PBMC_GENES = "CD3E, MS4A1, CD14, LYZ, NKG7"
BRAIN_GENES = "Snap25, Mbp, Gfap, Pcp4"

# JavaScript setter for any Shiny input: selectize, bootstrap-select (pickerInput), plain select, radio / checkbox groups,
# single checkbox, text / numeric / textarea; {"index": [i, ...]} selects options by position (empty option ignored)
JS_SET = r"""
([id, val]) => {
  const el = document.getElementById(id);
  if (!el) return 'missing';
  const $el = window.jQuery(el);
  const optionValues = () => Array.from(el.options || []).map(o => o.value).filter(v => v !== '');
  // {"index": [i]} picks by position; {"match": "regex"} picks the first option matching the regex (case-insensitive), else the first option
  const pick = (vals) => {
    if (val && typeof val === 'object' && !Array.isArray(val) && 'index' in val) return val.index.map(i => vals[i]).filter(v => v !== undefined);
    if (val && typeof val === 'object' && !Array.isArray(val) && 'match' in val) {
      // One regex or several in priority order: the first regex matching any option wins
      for (const m of [].concat(val.match)) { const re = new RegExp(m, 'i'); const hit = vals.find(v => re.test(v)); if (hit !== undefined) return hit; }
      return vals[0];
    }
    return val;
  };
  if (el.classList.contains('shiny-input-radiogroup') || el.classList.contains('shiny-input-checkboxgroup')) {
    const boxes = Array.from(el.querySelectorAll('input'));
    const wanted = [].concat(pick(boxes.map(b => b.value)));
    boxes.forEach(b => { b.checked = wanted.includes(b.value); });
    $el.find('input').first().trigger('change');
    return 'group:' + wanted.join('|');
  }
  if (el.tagName === 'SELECT') {
    if (el.selectize) {
      const sz = el.selectize;
      const vals = Object.keys(sz.options).filter(v => v !== '');
      const v = pick(vals);
      // Server-side selectize only holds the options already searched: add a missing value before selecting it
      [].concat(v).forEach(x => { if (typeof x === 'string' && !(x in sz.options)) { const o = {}; o[sz.settings.valueField] = x; o[sz.settings.labelField] = x; sz.addOption(o); } });
      sz.setValue(v, false);
      // Value actually held by the widget: an empty read-back means the selection did not take
      const got = [].concat(sz.getValue()).filter(x => x !== '');
      return got.length ? 'selectize:' + got.join('|') : 'notset:' + [].concat(v).join('|');
    }
    const v = pick(optionValues());
    if ($el.selectpicker && el.classList.contains('selectpicker')) { $el.selectpicker('val', v); $el.trigger('change'); return 'picker:' + [].concat(v).join('|'); }
    if (el.multiple) { const w = [].concat(v); Array.from(el.options).forEach(o => { o.selected = w.includes(o.value); }); }
    else { el.value = [].concat(v)[0]; }
    $el.trigger('change');
    return 'select:' + [].concat(v).join('|');
  }
  // sliderInput (ion.rangeSlider): a number, or [from, to] for a range slider
  const slider = $el.data('ionRangeSlider');
  if (slider) {
    const v = [].concat(val);
    slider.update(v.length > 1 ? {from: v[0], to: v[1]} : {from: v[0]});
    $el.trigger('change');
    return 'slider:' + v.join('|');
  }
  if (el.type === 'checkbox') { el.checked = !!val; $el.trigger('change'); return 'checkbox:' + el.checked; }
  $el.val(val); $el.trigger('change'); $el.trigger('input');
  return 'value:' + val;
}
"""

# State of the page after a step: errors to record, busy flag, open modal
JS_STATE = r"""
() => {
  const visible = (e) => !!(e.offsetWidth || e.offsetHeight || e.getClientRects().length);
  const notes = Array.from(document.querySelectorAll('#shiny-notification-panel .shiny-notification')).map(n => ({
    type: n.classList.contains('shiny-notification-error') ? 'error' : (n.classList.contains('shiny-notification-warning') ? 'warning' : 'message'),
    text: (n.querySelector('.shiny-notification-content-text') || n).innerText.trim()}));
  const pane = document.querySelector('.tab-pane.active') || document.body;
  const outErrors = Array.from(pane.querySelectorAll('.shiny-output-error')).filter(e => visible(e) && !e.classList.contains('shiny-output-error-validation'))
    .map(e => (e.id || '') + ': ' + e.innerText.trim().slice(0, 300));
  const modal = Array.from(document.querySelectorAll('.modal')).find(m => visible(m));
  return {
    busy: document.documentElement.classList.contains('shiny-busy'),
    disconnected: !!document.getElementById('shiny-disconnected-overlay'),
    notes: notes, outErrors: outErrors,
    modal: modal ? {text: modal.innerText.trim().slice(0, 500), hasFooterButton: !!modal.querySelector('.modal-footer button, .modal-footer a')} : null
  };
}
"""

# Visible action / download buttons with an id in the active tab pane (sweep)
JS_BUTTONS = r"""
() => {
  const visible = (e) => !!(e.offsetWidth || e.offsetHeight || e.getClientRects().length);
  const pane = document.querySelector('.tab-pane.active') || document.body;
  return Array.from(pane.querySelectorAll('button.action-button[id], a.action-button[id], a.shiny-download-link[id], button.shiny-download-link[id]'))
    .filter(e => visible(e) && !e.disabled && !e.closest('.modal'))
    .map(e => ({id: e.id, kind: e.classList.contains('shiny-download-link') ? 'download' : 'click', label: e.innerText.trim().slice(0, 60)}));
}
"""
# Buttons the sweep never clicks: they wipe the analysis, stop the app or only navigate
SKIP_SWEEP = re.compile(r"reset|delete|remove|clear|quit|stop|restart|exit|close|cancel|logout|load_workspace|restore|help|tutorial", re.I)


class Runner:
    def __init__(self, page, out_dir):
        self.page, self.out, self.results, self.module, self.tab = page, out_dir, [], "", ""
        self.tested, self.tabs = set(), {}  # ids already exercised by the scenarios; tabs visited per module (for the sweep)
        os.makedirs(os.path.join(out_dir, "screenshots"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "downloads"), exist_ok=True)

    def state(self):
        return self.page.evaluate(JS_STATE)

    def wait_idle(self, timeout_s=600, keep_modal=False):
        """Wait until Shiny is idle for 3 s: no busy class and no processing modal (a modal without footer buttons)"""
        t0, calm, idle_modal = time.time(), 0, 0
        while time.time() - t0 < timeout_s:
            s = self.state()
            if s["disconnected"]:
                return s
            processing = s["modal"] is not None and not s["modal"]["hasFooterButton"] and not keep_modal
            calm = calm + 1 if (not s["busy"] and not processing) else 0
            # A footer-less modal with R idle for 20 s is a form, not a processing message
            idle_modal = idle_modal + 1 if (processing and not s["busy"]) else 0
            if calm >= 3 or idle_modal >= 20:
                return s
            time.sleep(1)
        s = self.state()
        s["timeout"] = True
        return s

    def clear(self, keep_modal=False):
        """Remove notifications and close result / error modals so the next step starts clean"""
        self.page.evaluate("() => document.querySelectorAll('#shiny-notification-panel .shiny-notification').forEach(n => n.remove())")
        if not keep_modal:
            s = self.state()
            if s["modal"] is not None and s["modal"]["hasFooterButton"]:
                self.page.keyboard.press("Escape")
                self.page.evaluate("() => { const b = document.querySelector('.modal .modal-footer button'); if (b) b.click(); }")
                time.sleep(1)

    def record(self, step, status, t0, message=""):
        shot = ""
        if status in ("KO", "FATAL"):
            shot = f"screenshots/{len(self.results):03d}_{self.module}_{re.sub(r'[^A-Za-z0-9_]', '_', step)[:60]}.png"
            try:
                self.page.screenshot(path=os.path.join(self.out, shot), full_page=False)
            except Exception:
                shot = ""
        self.results.append({"module": self.module, "tab": self.tab, "step": step, "status": status, "start": t0, "end": time.time(),
                             "seconds": round(time.time() - t0, 1), "message": message[:1500], "screenshot": shot})
        print(f"[{status}] {self.module} / {self.tab} / {step} ({time.time() - t0:.1f}s) {message[:200]}", flush=True)
        with open(os.path.join(self.out, "ui_results.json"), "w") as f:
            json.dump(self.results, f, indent=1)

    def judge(self, step, t0, s, keep_modal=False):
        """Turn the page state after a step into OK / WARN / KO / FATAL"""
        if s.get("disconnected"):
            self.record(step, "FATAL", t0, "Shiny session disconnected (R crashed or the session ended): see app_log.txt")
            raise RuntimeError("disconnected")
        errors = [n["text"] for n in s["notes"] if n["type"] == "error"] + s["outErrors"]
        if s["modal"] and re.search(r"\berror\b", s["modal"]["text"], re.I):
            errors.append("Modal: " + s["modal"]["text"])
        warnings = [n["text"] for n in s["notes"] if n["type"] == "warning"]
        if s.get("timeout"):
            errors.append("Still busy after the time limit")
        if errors:
            self.record(step, "KO", t0, " || ".join(errors))
        elif warnings:
            self.record(step, "WARN", t0, " || ".join(warnings))
        else:
            self.record(step, "OK", t0, " | ".join(n["text"] for n in s["notes"])[:300])
        self.clear(keep_modal=keep_modal)

    # ---- actions ----
    def goto(self, module, tab):
        self.module, self.tab = module, tab
        self.tabs.setdefault(module, [])
        if tab not in self.tabs[module]:
            self.tabs[module].append(tab)
        t0 = time.time()
        found = self.page.evaluate("(t) => { const a = document.querySelector(`a[data-value='${t}']`); if (!a) return false; a.click(); return true; }", tab)
        if not found:
            self.record("open tab", "KO", t0, f"tab '{tab}' not found in the sidebar")
            return
        time.sleep(1.5)
        self.judge("open tab", t0, self.wait_idle(120))

    def set(self, input_id, value, quiet=False):
        t0 = time.time()
        res = self.page.evaluate(JS_SET, [input_id, value])
        time.sleep(0.8)
        if res == "missing":
            self.record(f"set {input_id}", "KO", t0, "input not found in the page (renamed or not rendered yet)")
        elif str(res).startswith("notset:"):
            self.record(f"set {input_id}", "KO", t0, f"value not accepted by the widget ({res[7:]}): option missing from the menu")
        else:
            s = self.wait_idle(120, keep_modal=True)
            if not quiet:
                self.judge(f"set {input_id} = {res}", t0, s, keep_modal=True)

    def upload(self, input_id, path, timeout_s=900, keep_modal=False):
        t0 = time.time()
        try:
            self.page.set_input_files(f"#{input_id}", path, timeout=60000)
        except Exception as e:
            self.record(f"upload {os.path.basename(path)} -> {input_id}", "KO", t0, str(e)[:500])
            return
        time.sleep(3)
        self.judge(f"upload {os.path.basename(path)} -> {input_id}", t0, self.wait_idle(timeout_s, keep_modal=keep_modal), keep_modal=keep_modal)

    def click(self, button_id, timeout_s=900, keep_modal=False):
        t0 = time.time()
        self.tested.add(button_id)
        loc = self.page.locator(f"#{button_id}")
        if loc.count() == 0:
            self.record(f"click {button_id}", "KO", t0, "button not found in the page")
            return
        if not loc.first.is_visible():
            self.record(f"click {button_id}", "SKIP", t0, "button hidden (conditional panel or collapsed box)")
            return
        if loc.first.is_disabled():
            self.record(f"click {button_id}", "KO", t0, "button disabled")
            return
        loc.first.scroll_into_view_if_needed()
        loc.first.click()
        time.sleep(2)
        self.judge(f"click {button_id}", t0, self.wait_idle(timeout_s, keep_modal=keep_modal), keep_modal=keep_modal)

    def download(self, link_id, timeout_s=600):
        t0 = time.time()
        self.tested.add(link_id)
        loc = self.page.locator(f"#{link_id}")
        if loc.count() == 0:
            self.record(f"download {link_id}", "KO", t0, "download button not found")
            return
        if not loc.first.is_visible():
            # Button inside an inactive tab of a tabsetPanel: open that tab first
            self.page.evaluate("""(id) => { const el = document.getElementById(id); const pane = el && el.closest('.tab-pane');
                if (pane && pane.id) { const a = document.querySelector(`a[href='#${pane.id}']`); if (a) a.click(); } }""", link_id)
            time.sleep(1.5)
        if not loc.first.is_visible():
            self.record(f"download {link_id}", "SKIP", t0, "download button hidden")
            return
        # Disabled download button: the result it exports does not exist yet (bind_download_state)
        if self.page.evaluate("(id) => { const el = document.getElementById(id); return !!el && (el.classList.contains('disabled') || el.hasAttribute('disabled')); }", link_id):
            self.record(f"download {link_id}", "SKIP", t0, "download disabled: no result yet")
            return
        try:
            loc.first.scroll_into_view_if_needed()
            with self.page.expect_download(timeout=timeout_s * 1000) as info:
                loc.first.click()
            dl = info.value
            dest = os.path.join(self.out, "downloads", f"{link_id}__{dl.suggested_filename}")
            dl.save_as(dest)
            size = os.path.getsize(dest)
            with open(dest, "rb") as f:
                head = f.read(200)
            if size == 0:
                self.record(f"download {link_id}", "KO", t0, "empty file")
            elif head.lstrip().lower().startswith(b"<!doctype html") or head.lstrip().lower().startswith(b"<html"):
                self.record(f"download {link_id}", "KO", t0, f"server returned an HTML page instead of a file ({dl.suggested_filename})")
            else:
                self.record(f"download {link_id}", "OK", t0, f"{dl.suggested_filename} ({size:,} bytes)")
        except PWTimeout:
            self.record(f"download {link_id}", "KO", t0, "no file received (handler error or nothing to download yet)")
        except Exception as e:
            self.record(f"download {link_id}", "KO", t0, str(e)[:500])
        self.clear()

    def sweep(self, module, timeout_s=300):
        """Click every visible action / download button of the module's tabs that the scenario did not exercise.
        Destructive or navigation-only buttons are skipped; two passes per tab catch buttons rendered by the first clicks."""
        for tab in list(self.tabs.get(module, [])):
            self.goto(module, tab)
            # Expand collapsed boxes so their buttons become visible
            self.page.evaluate("() => document.querySelectorAll('.tab-pane.active .box.collapsed-box [data-widget=collapse]').forEach(b => b.click())")
            time.sleep(1)
            for _ in range(2):
                found = self.page.evaluate(JS_BUTTONS)
                todo = [b for b in found if b["id"] not in self.tested]
                if not todo:
                    break
                for b in todo:
                    if SKIP_SWEEP.search(b["id"]) or SKIP_SWEEP.search(b["label"]):
                        self.tested.add(b["id"])
                        self.record(f"sweep {b['kind']} {b['id']}", "SKIP", time.time(), f"not clicked (destructive or navigation): '{b['label']}'")
                        continue
                    self.tab = tab
                    if b["kind"] == "download":
                        self.download(b["id"], timeout_s=timeout_s)
                    else:
                        self.click(b["id"], timeout_s=timeout_s)
                        # Close a form modal opened by the button (no footer button means wait_idle treated it as processing)
                        if self.state()["modal"] is not None:
                            self.page.keyboard.press("Escape")
                            self.page.evaluate("() => { const b = document.querySelector('.modal .modal-footer button, .modal button.close'); if (b) b.click(); if (window.jQuery) jQuery('.modal').modal('hide'); }")
                            time.sleep(1.5)
                    # A button may have switched tab: come back if the tab's sidebar link is no longer active
                    if not self.page.evaluate("(t) => { const a = document.querySelector(`a[data-value='${t}']`); return !!a && (a.classList.contains('active') || !!a.closest('li.active')); }", tab):
                        self.goto(module, tab)

    def snapshot(self, name):
        try:
            self.page.screenshot(path=os.path.join(self.out, "screenshots", f"view_{self.module}_{name}.png"), full_page=True)
        except Exception:
            pass


# ---------------------------------------------------------------------------------------------------------------------
# Module scenarios: real analysis order on the demonstration data

def single(r, demo):
    m = "single"
    r.goto(m, "load_dataset")
    r.set("species_choice", "human"); r.set("dataset_type", "snRNA")
    r.upload("file", f"{demo}/pbmc_1k_v3_filtered_feature_bc_matrix.h5")
    r.goto(m, "qc")
    # PBMC (single cells, not nuclei): the app defaults (mitochondrial < 5 %, 200-3500 genes) are meant for snRNA-seq and keep ~2 % of these cells
    r.set("nFeature_range", [200, 6000]); r.set("percent.mt_max", 20)
    for b in ["show_plots", "apply_qc", "normalize_data"]: r.click(b)
    r.goto(m, "dimensional_reduction")
    for b in ["scale_button", "run_elbow"]: r.click(b)
    r.goto(m, "clustering")
    for b in ["run_neighbors", "compute_3d_umap", "run_clustering"]: r.click(b)
    r.download("downloadUMAP"); r.snapshot("clustering")
    r.goto(m, "doublet_detection")
    r.click("run_doubletfinder", 1800); r.download("download_doublet_results"); r.click("remove_doublets")
    r.goto(m, "plots_genes_expressions")
    for i in ["gene_list_feature", "gene_list_dotplot", "gene_list_vln", "gene_list_ridge_plot"]: r.set(i, PBMC_GENES)
    r.set("gene_select_genes_analysis", ["CD3E", "MS4A1", "LYZ"])
    for b in ["show_feature", "show_dot", "show_vln", "show_ridge", "analyze_btn"]: r.click(b)
    for d in ["downloadFeaturePlot", "downloadDotPlot", "downloadVlnPlot", "download_ridge_plot", "download_genes_number_expresion"]: r.download(d)
    r.goto(m, "assigning_cell_type_identity")
    r.set("rename_single_cluster", "Renamed_test"); r.click("rename_single_cluster_button"); r.click("undo_cluster_single")
    r.set("singler_species_single", "human"); time.sleep(2)
    r.click("run_singler_single", 2400); r.click("apply_singler_labels_single")
    for d in ["download_singler_clusters_single", "download_singler_cells_single"]: r.download(d)
    r.click("update_colour"); r.snapshot("annotation")
    r.goto(m, "cluster_comparison")
    r.click("find_markers_global_single", 1200)
    r.set("cluster1_pairwise_single", {"index": [0]}); r.set("cluster2_pairwise_single", {"index": [0]})
    r.click("compare_markers_pairwise_single", 1200); r.click("generate_cluster_composition_single")
    for d in ["downloadUMAPCluster", "download_markers_global_single", "download_markers_pairwise_single", "download_cluster_composition_single"]: r.download(d)
    r.goto(m, "enrichment_single")
    r.set("enrichment_gene_sets_single", "human")
    r.set("enrichment_group1_single", {"index": [0]}); r.set("enrichment_group2_single", {"index": [1]})
    r.click("enrichment_run_single", 1800)
    for d in ["download_enrichment_ora_plot_single", "download_enrichment_ora_table_single", "download_enrichment_gsea_plot_single",
              "download_enrichment_gsea_table_single", "download_enrichment_params_single"]: r.download(d)
    r.goto(m, "exclusive_biomarkers")
    r.set("exclusive_target_cluster_single", {"index": [0]})
    r.click("find_exclusive_markers_single", 1200); r.download("download_exclusive_markers_single")
    r.set("exclusive_genes_to_plot_single", PBMC_GENES); r.click("generate_exclusive_plot_single"); r.download("download_exclusive_plot_single")
    r.goto(m, "deg_visualizations")
    r.set("venn_table_1_single", {"index": [0]}); r.set("venn_table_2_single", {"index": [1]})
    r.click("generate_venn_btn_single")
    for d in ["download_venn_diagram_single", "download_venn_gene_lists_single"]: r.download(d)
    r.set("volcano_deg_source_single", {"index": [0]}); r.click("generate_volcano_btn_single"); r.download("download_volcano_plot_single")
    r.goto(m, "heatmaps_dualexpression")
    r.set("gene_select_heatmap", PBMC_GENES); r.click("generateCustomHeatmap"); r.download("download_heatmap_single")
    r.set("gene_source_heatmap_single", "cluster_markers"); r.click("generateCustomHeatmap", 1200)
    r.download("download_heatmap_markers_single"); r.download("download_heatmap_single")
    r.goto(m, "coexpression_single")
    r.set("gene_text_coexpression_single", "CD3E, CD8A")
    r.set("feature1_select", "CD3E"); r.set("feature2_select", "CD8A"); r.click("generateScatter_single"); r.download("download_scatter_single")
    for unit in ["counts", "data"]:
        r.set("coexpr_threshold_unit_single", unit)
        r.set("coexpr_background_clusters_single", {"index": [0]})
        r.click("estimate_coexpr_thresholds_single"); r.click("analyze_coexpression_single")
    r.download("download_coexpression_table_single"); r.download("download_coexpression_plot_single")
    r.set("coexpr_de_cluster_single", {"index": [0]}); r.click("run_coexpr_de_single", 1200); r.download("download_coexpr_de_single")
    r.goto(m, "subset")
    r.set("select_ident_subset", {"index": [0, 1]}); r.click("apply_subset")
    r.set("gene_list_subset", "CD3E"); r.click("apply_gene_subset"); r.download("download_subset_seurat")
    r.goto(m, "session_single")
    for d in ["download_params_single", "download_report_single", "download_session_single"]: r.download(d, 1200)


# Human genes frequent in PBMC communication (MIF, MHC-II, integrins, chemokines): the list must hit significant interactions
PBMC_GENES = "MIF CD74 CXCR4 CD44 APP LGALS9 HLA-DRA HLA-DRB1 HLA-DPA1 CD99 ITGB2 ICAM1 ANXA1 FPR1 CCL5 CCR1 SELPLG SELL"


def cellchat(r, demo):
    m = "cellchat"
    r.goto(m, "load_data_cellchat")
    r.set("cellchat_import_source", "single"); r.click("import_session_object_cellchat")
    r.set("species_cellchat", "human")
    # Individual sources (all ticked), then single tiers, then the merged default used for the analysis
    r.set("db_mode_cellchat", "custom"); r.click("load_db_cellchat", 1200)
    r.set("db_mode_cellchat", "merged")
    r.set("db_filter_tiers_cellchat", ["high"]); r.click("load_db_cellchat")
    r.set("db_filter_tiers_cellchat", ["high", "medium", "low"]); r.click("load_db_cellchat")
    r.set("db_filter_tiers_cellchat", ["high", "medium"]); r.click("load_db_cellchat"); r.snapshot("database")
    # Genes of interest: a named list, replaced under the same name (one gene absent from the database), database query
    r.set("gene_list_name_cellchat", "PBMC test")
    r.set("gene_list_text_cellchat", PBMC_GENES); r.click("add_gene_list_cellchat")
    r.set("gene_list_text_cellchat", PBMC_GENES + " NOTAGENE1"); r.click("add_gene_list_cellchat")
    r.set("gene_list_selected_cellchat", "PBMC test"); r.click("run_gene_query_cellchat"); r.snapshot("gene_query")
    r.goto(m, "ligand_receptor_cellchat")
    r.set("cellchat_object_name", "pbmc_test")
    r.click("create_and_analyze_cellchat", 3600)
    # Bubble plot needs sources and targets (two first groups as senders, three first as receivers)
    r.set("sources_use_cellchat", {"index": [0, 1]}); r.set("targets_use_cellchat", {"index": [0, 1, 2]})
    r.click("generate_plot_cellchat", 600)
    for d in ["download_analyzed_objects", "download_bubble_plot", "download_info_table"]: r.download(d, 900)
    # Bubble plot restricted to the gene list, then back to all interactions
    r.set("bubble_gene_list_cellchat", "PBMC test"); r.click("generate_plot_cellchat", 600); r.snapshot("bubble_gene_list")
    r.download("download_info_table", 900)
    r.set("bubble_gene_list_cellchat", ""); r.click("generate_plot_cellchat", 600)
    r.goto(m, "circle_plot_cellchat")
    r.set("sender_groups", {"index": [0, 1]}); r.set("receiver_groups", {"index": [0, 1, 2]})
    r.click("generate_chord_plot", 900); r.download("download_chord_plot")
    # Chord plot restricted to the gene list (before the top-n selection), then back to all interactions
    r.set("chord_gene_list_cellchat", "PBMC test"); r.click("generate_chord_plot", 900); r.snapshot("chord_gene_list")
    r.set("chord_gene_list_cellchat", ""); r.click("generate_chord_plot", 900)
    r.set("split_chord_mode", "split")
    r.set("split_chord_senders", {"index": [0, 1]}); r.set("split_chord_receivers", {"index": [0, 1]})
    r.click("load_split_chord_lr"); r.click("generate_split_chord", 900); r.download("download_split_chord")


def trajectory(r, demo):
    m = "trajectory"
    r.goto(m, "load_trajectory")
    r.set("trajectory_import_source", "single"); r.click("import_session_object_trajectory", 900)
    r.goto(m, "trajectory_construction")
    r.click("constructGraph", 1800)
    r.set("root_cluster_select", {"index": [0]}); r.click("set_root_cell", 900)
    for d in ["download_trajectory_umap", "download_pseudotime_umap"]: r.download(d)
    r.goto(m, "differentialy_expressed_genes_trajectory")
    r.click("run_graph_test", 2400); r.download("download_graph_test")
    r.click("run_correlation_test", 1200); r.download("download_correlation")
    r.goto(m, "genes_visualization_trajectory")
    r.set("single_gene_umap_trajectory", "LYZ"); r.click("plot_single_gene_umap_trajectory"); r.download("download_single_gene_umap_trajectory")
    r.set("multi_gene_pseudotime_trajectory", ["LYZ", "CD3E"]); r.click("plot_multi_gene_pseudotime_trajectory")
    r.download("download_multi_gene_pseudotime_trajectory")


def multiple(r, demo):
    m = "multiple"
    r.goto(m, "load_datasets_merge")

    # Opening the loading modal resets every loaded dataset (new load cycle), so both files are loaded again before each method
    files = {1: "pbmc_1k_v3_filtered_feature_bc_matrix.h5", 2: "pbmc_1k_v2_filtered_feature_bc_matrix.h5"}

    def load_datasets():
        r.click("open_file_input_modal", 60, keep_modal=True)
        r.set("species_choice_merge", "human")
        # PBMC are single cells: the snRNA-oriented defaults (MT < 5 %, < 3500 genes) keep only ~2 % of the v3 cells
        r.set("min_features_merge", 200); r.set("max_features_merge", 6000); r.set("max_mt_percent_merge", 20)
        r.set("num_datasets", 2); time.sleep(2)
        for i, f in files.items():
            r.set(f"dataset_name{i}", ["PBMC_v3", "PBMC_v2"][i - 1]); r.set(f"dataset_type_merge{i}", "snRNA_merge")
            r.upload(f"merge{i}", f"{demo}/{f}", 1200, keep_modal=True)

    # Harmony last: its result is kept for the rest of the analysis
    for method in ["simple_merge", "integrate_fastmnn", "integrate_rpca", "integrate", "integrate_harmony"]:
        load_datasets()
        r.click(method, 3600)
        # A modal left open (failed method) would block every following click
        r.page.keyboard.press("Escape"); r.page.evaluate("() => { if (window.jQuery) jQuery('.modal').modal('hide'); }"); time.sleep(1)
    r.goto(m, "clustering_merge")
    for b in ["runScalePCA", "runFindNeighbors", "compute_3d_umap_merge", "runFindClusters"]: r.click(b, 1200)
    r.download("downloadUMAP_merge"); r.snapshot("clustering")
    r.goto(m, "plots_genes_expressions_merge")
    for i in ["gene_list_feature_merge", "gene_list_dot_merge", "gene_list_vln_merge", "gene_list_ridge_merge", "gene_list_genes_expression_merge", "gene_list_merge"]:
        r.set(i, PBMC_GENES, quiet=True)
    for b in ["show_feature_plot_merge", "runDotPlot", "runVlnPlot", "runRidgePlot", "analyze_btn_genes_expression_merge"]: r.click(b)
    for d in ["downloadFeaturePlotMerge", "downloadDotPlotMerge", "downloadVlnPlotMerge", "downloadRidgePlotMerge", "download_genes_number_expression_merge"]: r.download(d)
    r.goto(m, "assigning_cell_type_identity_merge")
    r.set("singler_species_merge", "human", quiet=True); time.sleep(2)
    r.click("run_singler_merge", 2400); r.click("apply_singler_labels_merge")
    for d in ["download_singler_clusters_merge", "download_singler_cells_merge"]: r.download(d)
    r.goto(m, "cluster_comparison_merge")
    for b in ["calculate_DE", "compare_clusters_button", "compare_datasets_button", "generate_cluster_table_multiple"]: r.click(b, 1200)
    for d in ["download_markers_single_cluster_merge", "download_markers_multiple_clusters_merge", "download_diff_dataset_cluster",
              "download_cluster_composition_multiple", "download_pie_chart_merge", "download_filtered_umap_plot"]: r.download(d)
    r.goto(m, "enrichment_merge")
    r.set("enrichment_gene_sets_merge", "human")
    r.set("enrichment_group1_merge", {"index": [0]}); r.set("enrichment_group2_merge", {"index": [1]})
    r.click("enrichment_run_merge", 1800)
    for d in ["download_enrichment_ora_table_merge", "download_enrichment_gsea_table_merge", "download_enrichment_params_merge"]: r.download(d)
    r.goto(m, "exclusive_biomarkers_merge")
    # Grouping column and target group are required by the observer; the cluster restriction is optional and left empty
    # Group by clusters (dataset vs dataset gives no exclusive marker: both PBMC samples hold the same cell types)
    r.set("exclusive_groupby_col_merge", {"match": ["^seurat_clusters$", "^ClusterIdents$", "cluster"]}); r.set("exclusive_target_groups_merge", {"index": [0]})
    r.click("find_exclusive_markers_merge", 1200); r.download("download_exclusive_markers_merge")
    r.set("exclusive_genes_to_plot_merge", PBMC_GENES); r.click("generate_exclusive_plot_merge"); r.download("download_exclusive_plot_merge")
    r.goto(m, "heatmaps_dualexpression_merge")
    r.set("gene_select_heatmap_multi", PBMC_GENES); r.click("generateCustomHeatmap_multi"); r.download("download_heatmap_multi")
    r.goto(m, "coexpression_merge")
    r.set("gene_text_coexpression_multiple", "CD3E, CD8A")
    for unit in ["counts", "data"]:
        r.set("coexpr_threshold_unit_multiple", unit)
        r.set("coexpr_background_clusters_multiple", {"index": [0]})
        r.click("estimate_coexpr_thresholds_multiple"); r.click("analyze_coexpression_multiple")
    r.download("download_coexpression_table_multiple")
    r.set("coexpr_de_cluster_multiple", {"index": [0]}); r.click("run_coexpr_de_multiple", 1200); r.download("download_coexpr_de_multiple")
    r.goto(m, "deg_visualizations_merge")
    r.set("venn_table_1", {"index": [0]}); r.set("venn_table_2", {"index": [1]}); r.click("generate_venn_btn")
    r.download("download_venn_diagram")
    r.set("volcano_deg_source_merge", {"index": [0]}); r.click("generate_volcano_btn_merge"); r.download("download_volcano_plot_merge")
    r.goto(m, "subset_merge")
    r.set("select_ident_subset_merge", {"index": [0, 1]}); r.click("apply_subset_merge"); r.download("download_subset_merge")
    r.goto(m, "session_merge")
    for d in ["download_params_multi", "download_report_multi", "download_session_merge"]: r.download(d, 1200)


def spatial(r, demo):
    m = "spatial"
    r.goto(m, "load_spatial_dataset")
    r.set("spatial_species_choice", "mouse")
    r.upload("spatial_file", f"{demo}/visium_mouse_brain.zip", 1800); r.snapshot("loaded")
    r.goto(m, "spatial_sections")
    for b in ["detect_sections", "apply_sections", "remove_debris_spots"]: r.click(b)
    r.goto(m, "spatial_normalisation")
    for b in ["spatial_qc_plots", "apply_spatial_qc", "normalize_spatial_data"]: r.click(b, 1800)
    r.goto(m, "spatial_clustering")
    for b in ["spatial_scale_pca", "spatial_run_neighbors_umap", "spatial_find_clusters"]: r.click(b, 1800)
    r.download("save_spatial_object", 1200); r.snapshot("clustering")
    r.goto(m, "spatial_interactive_visualization")
    for b in ["viewer_zoom_in", "viewer_pan_left", "viewer_pan_up", "viewer_pan_down", "viewer_pan_right", "viewer_zoom_out", "viewer_reset"]: r.click(b)
    r.set("viewer_color_by", "gene"); r.set("viewer_gene", "Mbp"); r.download("download_viewer_plot"); r.snapshot("viewer")
    r.goto(m, "spatial_gene_expression_visualisation")
    for i in ["spatial_genes_text", "spatial_violin_genes_text", "spatial_dotplot_genes_text"]: r.set(i, BRAIN_GENES)
    for b in ["show_spatial_features", "show_spatial_violin", "show_spatial_dotplot"]: r.click(b)
    for d in ["download_spatial_features", "download_spatial_violin", "download_spatial_dotplot"]: r.download(d)
    r.goto(m, "spatial_cluster_annotation")
    r.click("apply_all_spatial_names")
    for d in ["download_spatial_umap", "download_spatial_tissue"]: r.download(d)
    r.goto(m, "spatial_colocalization")
    r.set("coloc_genes", "Snap25, Mbp"); r.set("coloc_n_permutations", 50)
    r.click("run_coloc_analysis", 1800)
    for d in ["download_coloc_results", "download_coloc_map", "download_coloc_bar"]: r.download(d)
    r.goto(m, "spatial_svg")
    r.set("svg_n_features", 500); r.click("run_svg", 1800)
    for d in ["download_svg_table", "download_svg_plot"]: r.download(d)
    r.goto(m, "spatial_marker_analysis")
    r.set("spatial_marker_cluster", {"index": [0]}); r.click("run_spatial_markers", 1200); r.download("download_spatial_markers")


MODULES = {"single": single, "cellchat": cellchat, "trajectory": trajectory, "multiple": multiple, "spatial": spatial}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://cellhub_app:3838")
    ap.add_argument("--demo", default="/demo")
    ap.add_argument("--out", default="/out")
    ap.add_argument("--modules", default="single,cellchat,trajectory,multiple,spatial")
    ap.add_argument("--no-sweep", dest="no_sweep", action="store_true", help="only the scenarios, no sweep of the remaining buttons")
    a = ap.parse_args()
    with sync_playwright() as p:
        browser = p.chromium.launch(args=["--disable-dev-shm-usage"])
        ctx = browser.new_context(viewport={"width": 1600, "height": 1000}, accept_downloads=True)
        page = ctx.new_page()
        page.set_default_timeout(60000)
        # Browser console errors (JavaScript) kept for the report
        console = []
        page.on("console", lambda msg: console.append(f"{msg.type}: {msg.text}") if msg.type == "error" else None)
        r = Runner(page, a.out)
        single_done = False
        for name in a.modules.split(","):
            # Each module starts from a fresh page (new Shiny session), except CellChat and trajectory after the single module:
            # they import the single-dataset object, which only exists in the session where it was built
            try:
                if not (name in ("cellchat", "trajectory") and single_done):
                    page.goto(a.url, wait_until="networkidle", timeout=180000)
                r.module, r.tab = name, "start"
                r.wait_idle(120)
                MODULES[name](r, a.demo)
                if name == "single":
                    single_done = True
                if not a.no_sweep:
                    r.module = name
                    r.sweep(name)
            except Exception as e:
                r.module = name
                if str(e) != "disconnected":
                    r.record("module aborted", "FATAL", time.time(), traceback.format_exc()[-1500:])
        with open(os.path.join(a.out, "browser_console.txt"), "w") as f:
            f.write("\n".join(console))
        browser.close()


if __name__ == "__main__":
    main()
