Primer Design App

A lightweight web app for designing PCR primers (with an miRNA focus). Built with a simple Flask + Jinja stack so you can run it locally or deploy it anywhere a WSGI app runs.

✨ Features
Upload or paste sequences and get suggested primer pairs
Basic primer sanity checks (length, GC%, Tm, homopolymers)
Optional lookup against a bundled mirna_database.txt
Clean, minimal UI (Jinja templates + static assets)
One-file backend (app.py) that’s easy to read and extend

🧰 Tech Stack
Backend: Python (Flask)
Frontend: HTML/Jinja templates + vanilla JS/CSS
Data: Plain-text mirna_database.txt (tab/line-based)

🚀 Quickstart
1) Clone and install
<pre> 
  git clone https://github.com/Rusklass/primer_design_app.git
  cd primer_design_app
  
  # create a virtual environment
  python -m venv .venv
  # Windows: 
  .venv\Scripts\activate
  # macOS/Linux: 
  source .venv/bin/activate
  
  # install dependencies
  pip install -r requirements.txt

</pre>
2) Run the dev server
<pre> 
  # set Flask env (Windows PowerShell)
  $env:FLASK_APP="app.py"; $env:FLASK_ENV="development"
  
  # macOS/Linux
  export FLASK_APP=app.py
  export FLASK_ENV=development
  
  flask run 
</pre>
Open http://127.0.0.1:5000 in your browser.
  If flask isn’t available as a CLI, you can also do:
<pre> 
  python app.py
</pre>

📁 Project structure
<pre>
  primer_design_app/
├─ app.py                 # Flask application (routes + primer logic/glue)
├─ requirements.txt       # Python dependencies
├─ mirna_database.txt     # Optional reference data (miRNA entries)
├─ templates/             # Jinja2 HTML templates
└─ static/                # CSS/JS/assets
</pre>

🧪 Usage
Start the server.
Paste a target sequence or select a record (if the UI exposes the miRNA DB).
Adjust design parameters (if available: target length, GC%, Tm, etc.).
Submit to receive primer candidates and basic stats.
Exactly which inputs/outputs are shown depends on the current templates; the backend is kept minimal so you can easily adapt the form fields.

🧬 Primer design rules (baseline)
If you’re extending the logic, these common defaults are a good starting point:
Length: 18–24 nt
GC content: 40–60%
Melting temperature (Tm): 55–65 °C (pair within ~2 °C)
Avoid long homopolymers (e.g., >4 identical bases)
Check for hairpins/self-dimers/cross-dimers (simple heuristics or introduce a library)
You can keep this simple (string rules) or wire in a library for more rigorous checks.

⚙️ Configuration
Create an .env (optional) if you want to override defaults:
<pre>
  # .env
  FLASK_ENV=development
  SECRET_KEY=change-me
</pre>
Then load it in app.py (if not already) using python-dotenv:
<pre>
  from dotenv import load_dotenv
  load_dotenv()
</pre>
📦 Dependencies
Install via requirements.txt. If you add libraries (e.g., biopython, python-dotenv), remember to:
<pre>
  pip install <libname>
  pip freeze > requirements.txt
</pre>
🧯 Troubleshooting
Port already in use: set FLASK_RUN_PORT=5050 (or run flask run -p 5050).
Can’t import Flask: pip install -r requirements.txt in the active virtual env.
Templates not updating: ensure FLASK_ENV=development or restart the server.

🧩 Extending the app
Add more parameters to the form in templates/.
Encapsulate primer heuristics into helper functions/classes in app.py or a new services/ module.
Add export (CSV/FASTA) and basic report pages.
Integrate BLAST or off-target checks (server-side queue recommended).

🏗️ Production deployment
Gunicorn + gevent (Linux):
<pre>
  pip install gunicorn gevent
  gunicorn -k gevent -w 2 -b 0.0.0.0:8000 app:app
</pre>
Put Nginx in front for TLS and static file caching.
Set FLASK_ENV=production, SECRET_KEY to a strong random value.

✅ License
Copyright © 2025 Ruslan Klassen. All rights reserved.

🙌 Acknowledgements
Primer design heuristics are based on standard PCR guidelines.
miRNA dataset credit: miRBase [RELEASE 22.1] https://www.mirbase.org/.
