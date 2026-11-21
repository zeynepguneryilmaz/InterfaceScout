# Windows installation and launch

InterfaceScout v1.1.0-rc1 uses Python **3.10–3.12**. Extract the **Windows** Release ZIP completely into a normal directory (not a temporary ZIP preview).

## First-time setup

Double-click `run_local.bat` or run it in Command Prompt:

```bat
run_local.bat
```

The script creates `backend\.venv` if it does not already exist, checks the exact pinned dependencies and installs only if they are missing or incompatible. If installation fails, the command window stays open and explains where to look. For subsequent runs, use `start.bat`, which does not contact PyPI.

## CERTIFICATE_VERIFY_FAILED or self-signed certificate

The secure HTTPS connection used by pip cannot verify a certificate from the current network. Corporate/university HTTPS interception and custom trust-store configuration are common causes, but the specific origin must be confirmed locally.

**Do not use `--trusted-host`, `--cert ""`, `PYTHONHTTPSVERIFY=0`, or disable certificate verification.**

1. First check whether the same issue occurs on a different trusted network (e.g., a personal hotspot), if your institution permits that test.
2. If your organization manages the network, request its **approved CA certificate bundle in PEM format** from IT. The PEM bundle must include all necessary trusted roots.
3. From the same Command Prompt, set `PIP_CERT` and rerun the launcher:

```bat
set "PIP_CERT=C:\path\to\verified-CA-bundle.pem"
run_local.bat
```

Replace that example with the actual verified certificate path. `PIP_CERT` is local to this command window; it is not bundled with the application.

If the machine cannot access PyPI or institutional CA support is unavailable, use an IT-approved **offline wheelhouse** matching Windows, your Python version and architecture. Example commands:

```bat
backend\.venv\Scripts\python.exe -m pip install --no-index --find-links "C:\approved_wheels" -r backend\requirements.txt
start.bat
```

Use `run_local.bat` instead of `start.bat` if the environment has not been fully initialized. If a previous attempt failed, you generally do **not** need to delete the newly created virtual environment: the installer can retry.

## Troubleshooting

- `fastapi==0.115.12` is a pinned backend dependency. A download-time SSL error does not indicate a problem with FastAPI itself or InterfaceScout's chemistry calculations.
- To check the installed packages without downloading anything:
  `backend\.venv\Scripts\python.exe backend\dependency_check.py`
- `start.bat` logs startup errors to `backend\startup_log.txt`.
- The local interface opens at http://localhost:8000.

See [examples](../examples/README.md) for verified protein demonstrations.
