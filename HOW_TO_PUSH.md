# Pushing this folder to the existing repository

Run these on your own machine (authenticate with `gh auth login` or an SSH key; never paste a token into a chat or a file).

```bash
git clone https://github.com/vishal1601-2005/Entropy-walls-LVS.git
cd Entropy-walls-LVS
git tag v1-original && git push origin v1-original      # keeps the version the first manuscript used
git checkout -b revision-2026-10
rsync -a --exclude HOW_TO_PUSH.md --exclude .git /path/to/this/folder/ ./
git add -A
git status                                              # check what is going in
git commit -m "Revised manuscript: M_* Affleck-Dine scan, cross-model table, explicit amplitude check, LVS decay widths"
git push -u origin revision-2026-10
```
Open a pull request from `revision-2026-10` and merge when you are happy. Update the arXiv/Zenodo record or tag (e.g. `v2`) at submission so the repository history identifies the version used in the paper.
