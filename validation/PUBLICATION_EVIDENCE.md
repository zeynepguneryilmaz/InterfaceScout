# InterfaceScout v1.1.0 — publication evidence

Scientific model: `2.5.0-dual-accessibility`. Exact application source: `636bc5abfd72798965465ccb7bf79deb0ed2f951`. [Publication release](https://github.com/zeynepguneryilmaz/InterfaceScout/releases/tag/v1.1.0-npbe-platform-20261125). The eleven canonical maps use the same chemistry and geometric source files as the frozen original production revision `4c0943afbd3e8506c32925c82fcf4b33efc95209`.

| Evidence | Scope | Run |
|---|---|---|
| TV48 canonical verification | 48 structures, 528 maps, 15,910 regions | [data](FINAL_V110_RESULTS.md) |
| ISV76 canonical maps | 76 structures, 836 maps, 64,592 regions | [45385665397](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45385665397) |
| Browser atlas | 76 structures, original LPBE annotation where available | [45399987636](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45399987636) |
| Independent FreeSASA | 76 structures | [45445684505](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45445684505) |
| Final-model numerical sensitivity | 11,616 records | [45636665184](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45636665184) (sensitivity jobs; paired PB rerun separately) |
| Region-definition robustness | 5,376 records | [45663104217](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45663104217) |
| Paired LPBE/NPBE | 76 structures, 35,140 residue pairs | [45653251591](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45653251591) |
| 8C50 numerical grid diagnostic | 12 solver/grid configurations | [45657476304](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45657476304) |
| Paired APBS solver performance | 152 measurements | [45663636698](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45663636698) |
| Real solver regression | NPBE/LPBE/paired; canonical map invariance | [45657475939](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45657475939) |
| Platform tests and packages | Windows x64, macOS Intel, Linux x64 | [45661100267](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45661100267), [45664039001](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45664039001) |

The [persistent research data packages](../datasets/README.md) preserve source-run outputs and their file-level SHA-256 inventories. Actions links above identify provenance; the browser atlas remains an optional, retention-limited artifact. The archive manifest records source-run conclusions, including partial or failed external-server runs. 1GCT canonical maps use the original protein input; the paired PB calculations explicitly omit two N-only fragments during PQR preparation. The 8C50 numerical warnings are retained. Literature-associated cases provide application context for the descriptive maps.
