# Windows installation and launch

Supported Python versions: **3.10–3.12**. Open a Command Prompt in the repository root.

## First-time setup
```bat
run_local.bat
```

This creates `backend\.venv` and installs pinned dependencies. Later, launch with:
```bat
start.bat
```

Open http://localhost:8000 if the browser does not open automatically.

## If pip reports CERTIFICATE_VERIFY_FAILED

This is a trust-chain problem between your Python/pip installation and the package index, often caused by a corporate/university HTTPS inspection certificate. **Do not use --trusted-host or disable certificate verification.**

1. If a university or company manages your computer, ask IT for its **root CA certificate PEM file** and verify its provenance.
2. Use that certificate for pip without disabling verification:
   ```bat
   set "PIP_CERT=C:\path\to\trusted-university-ca-bundle.pem"
   run_local.bat
   ```
3. If `PIP_CERT` is insufficient, ensure the PEM bundle includes both the corporate root and the public trust roots appropriate to your system. Alternatively, request an IT-approved offline wheelhouse and install with `python -m pip install --no-index --find-links C:\wheelhouse -r backend\requirements.txt` in the activated virtual environment.
4. If already installed, use `start.bat`; it does **not** install packages.

Do not copy example export files into `backend/`. See [examples](../examples/README.md).
