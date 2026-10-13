# Bundled R package sources

These source archives are installed offline into a task-local R library. They
replace image-build-time `install_github()` calls and never write to a user R
library.

| Archive | Upstream | Commit |
|---|---|---|
| `presto-3c97180d.tar.gz` | `https://github.com/immunogenomics/presto` | `3c97180d90d330524b5b8353a342ab2f4b990d21` |
| `scmageck-e24aa0f3.tar.gz` | `https://github.com/weili-lab/scMAGeCK` | `e24aa0f3c9734c26ffbc200107adf0e199b26208` |

The archives were produced with `git archive` from the exact commits above;
their upstream license files are included inside each archive.

SHA-256:

```text
0293826e79eb905c16bfeeb5142b90f17b60b8c85f1441fa8be66700e6c75193  presto-3c97180d.tar.gz
f2179fc78991f780536c715b1f0dfdc21f48aae299767e17efcdd5b86b12acef  scmageck-e24aa0f3.tar.gz
```
