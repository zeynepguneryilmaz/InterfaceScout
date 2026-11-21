# Third-party software notices

InterfaceScout's original source code is licensed under the [MIT License](LICENSE). That license does **not** replace the licensing terms of redistributed third-party software.

## Bundled browser libraries

| Component | Included file | Version | License | Original project |
| --- | --- | --- | --- | --- |
| 3Dmol.js | `frontend/vendor/3Dmol-2.1.0.min.js` | 2.1.0 | BSD-3-Clause, with embedded-component notices | https://github.com/3dmol/3Dmol.js |
| SheetJS Community Edition | `frontend/vendor/xlsx-0.20.3.full.min.js` | 0.20.3 | Apache-2.0 | https://sheetjs.com/ |

**3Dmol.js**: Copyright (c) 2014, University of Pittsburgh and contributors. 3Dmol.js incorporates components from GLmol, Three.js and jQuery, with additional attribution and licensing terms reproduced in [the upstream license and incorporated-component notices](licenses/3Dmol.js-LICENSE.txt). The bundled minified file refers to a license-information file; this repository supplies the corresponding upstream notices here.

**SheetJS Community Edition**: Copyright (C) 2012–present SheetJS LLC. Licensed under the Apache License, Version 2.0. Attribution required by the upstream project:

```text
SheetJS Community Edition -- https://sheetjs.com/
Copyright (C) 2012-present   SheetJS LLC
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
      http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```

The [full Apache-2.0 license text](licenses/Apache-2.0.txt) is provided in this repository. These upstream licenses and notices must be retained when redistributing their respective components.

## Python dependencies

InterfaceScout depends on independently licensed Python packages specified in `backend/requirements.txt`, including FastAPI, Uvicorn, Pydantic, NumPy, Biopython, ReportLab, OpenPyXL, and Pillow. Those packages are installed separately by the user or provided as dependency wheels in applicable distributions; their own license notices govern those packages. Their versions are pinned in the requirements file.

## Trademarks and endorsement

Project names and logos are used solely to identify the upstream dependencies. Their inclusion does not imply sponsorship or endorsement of InterfaceScout by the upstream authors or institutions.

For full license conditions, consult the relevant license texts and upstream notices. This file is an attribution record, not a substitute for the complete license terms.
