# macOS installation and launch

Python **3.10–3.12** is required.

```bash
bash run_local.sh
```

For subsequent runs:
```bash
bash start.command
```

The local interface is at http://localhost:8000. A downloaded `.command` file may need permission: `chmod +x start.command`.

If pip reports a certificate-chain error, install your institution's verified CA into the trust chain and configure `PIP_CERT=/absolute/path/to/ca-bundle.pem` before `bash run_local.sh`. Do not disable TLS verification. For offline installation, use an approved local wheelhouse. Once installed, `bash start.command` reuses the existing environment.

See [examples](../examples/README.md).
