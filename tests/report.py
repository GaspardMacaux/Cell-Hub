"""Build the HTML report of a Cell-Hub test run from /out: ui_results.json, db_results.csv, mem.csv, app_log.txt"""
import csv, html, json, os, re, sys
out = sys.argv[1] if len(sys.argv) > 1 else "/out"
esc = html.escape
ui = json.load(open(os.path.join(out, "ui_results.json"))) if os.path.exists(os.path.join(out, "ui_results.json")) else []
db = list(csv.DictReader(open(os.path.join(out, "db_results.csv")))) if os.path.exists(os.path.join(out, "db_results.csv")) else []
fn = list(csv.DictReader(open(os.path.join(out, "functional_results.csv")))) if os.path.exists(os.path.join(out, "functional_results.csv")) else []
mem = []
if os.path.exists(os.path.join(out, "mem.csv")):
    for row in csv.reader(open(os.path.join(out, "mem.csv"))):
        try:
            mem.append((float(row[0]), float(row[1])))
        except (ValueError, IndexError):
            pass
log = open(os.path.join(out, "app_log.txt"), errors="replace").read().splitlines() if os.path.exists(os.path.join(out, "app_log.txt")) else []

def peak_gb(t0, t1):
    vals = [m for t, m in mem if t0 - 2 <= t <= t1 + 2]
    return f"{max(vals):.1f}" if vals else ""

color = {"OK": "#d4edda", "WARN": "#fff3cd", "KO": "#f8d7da", "FATAL": "#f5c6cb", "SKIP": "#e2e3e5", "INFO": "#d1ecf1"}
parts = ["<!doctype html><meta charset='utf-8'><title>Cell-Hub test report</title><style>body{font-family:sans-serif;margin:20px}"
         "table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #ccc;padding:4px 6px;vertical-align:top}th{background:#eee}"
         ".msg{max-width:700px;white-space:pre-wrap}</style><h1>Cell-Hub test report</h1>"]
# Summary per module
mods = {}
for r in ui:
    mods.setdefault(r["module"], {}).setdefault(r["status"], 0)
    mods[r["module"]][r["status"]] += 1
parts.append("<h2>Interface: summary per module</h2><table><tr><th>Module</th>" + "".join(f"<th>{s}</th>" for s in color) + "<th>Duration (min)</th><th>Peak RAM (GB)</th></tr>")
for mname, counts in mods.items():
    rs = [r for r in ui if r["module"] == mname]
    t0, t1 = min(r["start"] for r in rs), max(r["end"] for r in rs)
    parts.append(f"<tr><td>{esc(mname)}</td>" + "".join(f"<td>{counts.get(s, '')}</td>" for s in color) + f"<td>{(t1 - t0) / 60:.1f}</td><td>{peak_gb(t0, t1)}</td></tr>")
parts.append("</table>")
# Problems first, then every step
problems = [r for r in ui if r["status"] in ("KO", "FATAL", "WARN")]
def step_table(rows, title):
    parts.append(f"<h2>{title} ({len(rows)})</h2><table><tr><th>#</th><th>Module</th><th>Tab</th><th>Step</th><th>Status</th><th>s</th><th>RAM GB</th><th>Message</th><th>Screenshot</th></tr>")
    for i, r in enumerate(rows):
        shot = f"<a href='{esc(r['screenshot'])}'>png</a>" if r.get("screenshot") else ""
        parts.append(f"<tr style='background:{color.get(r['status'], '#fff')}'><td>{i}</td><td>{esc(r['module'])}</td><td>{esc(r['tab'])}</td><td>{esc(r['step'])}</td>"
                     f"<td>{r['status']}</td><td>{r['seconds']}</td><td>{peak_gb(r['start'], r['end'])}</td><td class='msg'>{esc(r['message'])}</td><td>{shot}</td></tr>")
    parts.append("</table>")
step_table(problems, "Interface: problems (KO / FATAL / WARN)")
# Functional test: app functions vs independent reference, in run order
if fn:
    n = {k: sum(1 for r in fn if r["status"] == k) for k in ("OK", "WARN", "KO", "SKIP")}
    parts.append(f"<h2>Functional test (app vs independent reference): OK {n['OK']} | WARN {n['WARN']} | KO {n['KO']} | SKIP {n['SKIP']}</h2>"
                 "<table><tr><th>Section</th><th>Check</th><th>Status</th><th>App value</th><th>Expected / reference</th><th>Detail</th></tr>")
    for r in fn:
        parts.append(f"<tr style='background:{color.get(r['status'], '#fff')}'><td>{esc(r['section'])}</td><td>{esc(r['check'])}</td><td>{r['status']}</td>"
                     f"<td>{esc(r['value'])}</td><td>{esc(r['expected'])}</td><td class='msg'>{esc(r['detail'])}</td></tr>")
    parts.append("</table>")
# Database checks
if db:
    parts.append("<h2>Database checks</h2><table><tr><th>Species</th><th>Check</th><th>Status</th><th>n</th><th>Detail</th></tr>")
    order = {"KO": 0, "WARN": 1, "INFO": 2, "OK": 3}
    for r in sorted(db, key=lambda r: (order.get(r["status"], 4), r["species"])):
        parts.append(f"<tr style='background:{color.get(r['status'], '#fff')}'><td>{esc(r['species'])}</td><td>{esc(r['check'])}</td><td>{r['status']}</td><td>{esc(r['n'])}</td><td class='msg'>{esc(r['detail'])}</td></tr>")
    parts.append("</table>")
# R console errors and warnings of the app container
err_lines = [l for l in log if re.search(r"error|erreur|warning|avis|fatal", l, re.I)]
parts.append(f"<h2>App console: error / warning lines ({len(err_lines)})</h2><pre style='font-size:12px;white-space:pre-wrap'>{esc(chr(10).join(err_lines[-400:]))}</pre>")
step_table(ui, "Interface: every step")
open(os.path.join(out, "report.html"), "w").write("\n".join(parts))
ko = sum(1 for r in ui if r["status"] in ("KO", "FATAL"))
print(f"Report: {os.path.join(out, 'report.html')} | interface steps {len(ui)}, KO/FATAL {ko} | database rows {len(db)}, KO {sum(1 for r in db if r['status'] == 'KO')} | functional checks {len(fn)}, KO {sum(1 for r in fn if r['status'] == 'KO')}")
