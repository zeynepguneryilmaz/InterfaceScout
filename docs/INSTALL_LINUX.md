# Linux installation and launch

Python **3.10–3.12** is required.

```bash
bash run_local.sh
```

For later runs:
```bash
bash start.sh
```

Open http://localhost:8000 if the browser does not open.

If pip fails with a self-signed certificate chain, consult your administrator for a verified CA bundle and run `PIP_CERT=/absolute/path/to/ca-bundle.pem bash run_local.sh`. Do not bypass TLS validation. For offline systems, install from a trusted wheelhouse into `backend/.venv`. Once the dependencies have been installed, `bash start.sh` does not invoke pip.

See [examples](../examples/README.md).

## Linux desktop icon

After extracting the Linux ZIP, run `bash "Create Desktop Shortcut.sh"`. This creates an `InterfaceScout.desktop` launcher with the new icon on your desktop and in your local applications directory. Some desktop environments require selecting **Allow Launching** once. Keep the extracted folder in place because the shortcut points to it.

