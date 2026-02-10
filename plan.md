# CRAN Readiness Plan for `mrsensemakr`

## Current State
The package is functionally sound with good test coverage and a working core. However, several issues need to be resolved before a CRAN submission will succeed. Infrastructure items are modeled after what's already in place for `sensemakr`.

---

## Critical Issues (will cause `R CMD check` ERRORs or WARNINGs)

### 1. Fix DESCRIPTION file
**File:** `DESCRIPTION`

- **Replace `Author`/`Maintainer` with `Authors@R`**: CRAN now requires the `Authors@R` field using `person()`. The current `Author: Carlos Cinelli,` (with trailing comma) and separate `Maintainer:` field are deprecated. (sensemakr uses `Authors@R` with multiple authors/roles.)
- **Remove trailing space from `Title`**: `"Sensitivity Analysis Tools for Mendelian Randomization "` has a trailing space — `R CMD check` will flag this.
- **Improve `Description`**: The current text starts with "The R package mrsensemakr" — CRAN policy says the Description field should not start with the package name or "This package". Rewrite to begin with the purpose.
- **Add `URL` and `BugReports`**: CRAN strongly recommends these. Point to the GitHub repo (`https://github.com/carloscinelli/mrsensemakr`). sensemakr uses both.
- **Update R version dependency**: `Depends: R (>= 2.10)` is extremely old. sensemakr uses `R (>= 3.1.0)`. Update to at least `R (>= 3.5.0)`.

### 2. Add missing roxygen2 documentation for exported S3 methods
**Files:** `R/mrsensemakr.R` (line 221), `R/mr-plots.R` (line 2)

- `print.mr_sensemakr` has only `##'@export` — needs `@param`, `@return`, `@description`.
- `plot.mr_sensemakr` has only `##'@export` — needs `@param`, `@return`, `@description`.
- The main `mr_sensemakr()` function's roxygen block is missing `@return` (describes what the function returns). This will generate a WARNING.

### 3. Fix `\value` (return value) documentation for all exported functions
`R CMD check` now warns on missing `\value` sections in `.Rd` files. Add `@return` tags to:
- `mr_sensemakr()`
- `print.mr_sensemakr()`
- `plot.mr_sensemakr()`

### 4. Improve data documentation
**File:** `R/data-documentation.R`

- `sim_data.Rd` needs a `@description` that is more than just repeating the title ("Simulated Data" for both title and description).
- The documentation is very sparse — just "A data frame with 200,000 observations and 27 variables." CRAN reviewers may request descriptions of the variables (at least the key ones: `out.trait`, `exp.trait`, `prs`, `age`, `sex`, `alcohol`, `smoking`, `pc1`-`pc20`). sensemakr documents all variables in its datasets.

### 5. Fix vignette: `T`/`F` instead of `TRUE`/`FALSE`
**File:** `vignettes/simulations.Rmd` (lines 90, 106, 112, 197)

Uses `F` and `T` shorthand (e.g., `R1.a = F`, `OUTLIERtest = T`). Even though the vignette is `eval=FALSE`, CRAN checks can flag this. Use `TRUE`/`FALSE` everywhere.

### 6. Fix test file: `rm(list = ls())` and `file.remove("Rplots.pdf")`
**File:** `tests/testthat/test-mrsensemakr.R`

- `rm(list = ls())` at the start of each test block is unnecessary in testthat (each test has its own environment) and is bad practice.
- `file.remove("Rplots.pdf")` at line 240 will fail if the file doesn't exist (e.g., on systems without a display). Replace with `unlink("Rplots.pdf")` which silently ignores missing files.

---

## CI/CD & Infrastructure (matching sensemakr)

### 7. Add GitHub Actions for R CMD check
**Create:** `.github/workflows/R-CMD-check.yaml`

sensemakr uses `r-lib/actions` with a matrix of 8 OS/R-version combos:
- macOS-latest (release)
- windows-latest (release, 4.1)
- ubuntu-latest (devel, release, oldrel-1, oldrel-2, oldrel-3, oldrel-4)

Create the same workflow for mrsensemakr. Key settings from sensemakr:
- Trigger on push to `[main, master]` and pull requests against `[main, master]`
- `fail-fast: false`
- Build args: `'--no-manual', '--compact-vignettes=gs+qpdf'`
- Uses: `r-lib/actions/setup-pandoc`, `r-lib/actions/setup-r`, `r-lib/actions/setup-r-dependencies`, `r-lib/actions/check-r-package`

### 8. Add GitHub Actions for test coverage (Codecov)
**Create:** `.github/workflows/test-coverage.yaml`

sensemakr runs coverage on ubuntu-latest, generates Cobertura XML, and uploads via `codecov/codecov-action@v4`. Replicate for mrsensemakr.

### 9. Update README badges
**File:** `README.Rmd` (and regenerate `README.md`)

Replace old Travis CI and AppVeyor badges with the same badge set sensemakr uses:
- **CRAN version badge**: `https://www.r-pkg.org/badges/version/mrsensemakr` linking to CRAN
- **CRAN downloads badge**: `https://cranlogs.r-pkg.org/badges/mrsensemakr`
- **R-CMD-check badge**: GitHub Actions status badge
- **Codecov badge**: `https://codecov.io/gh/carloscinelli/mrsensemakr/branch/master/graph/badge.svg`

Remove Travis and AppVeyor badges entirely.

### 10. Remove legacy CI config files
**Delete:** `.travis.yml`, `appveyor.yml`, `codecov.yml`

These are obsolete now that GitHub Actions handles CI and coverage. sensemakr has already deprecated these (they're only in `.Rbuildignore` as historical references).

---

## Important Issues (NOTEs / best practice)

### 11. Create `NEWS.md`
CRAN reviewers expect a changelog. Create a `NEWS.md` with at least an entry for the current version. sensemakr maintains a NEWS.md.

### 12. Create `inst/CITATION`
The package accompanies a published paper (Cinelli et al.). Create a proper citation file so `citation("mrsensemakr")` returns the correct reference. sensemakr has one.

### 13. Create `cran-comments.md`
**Create:** `cran-comments.md`

sensemakr maintains a `cran-comments.md` documenting test environments and submission notes for CRAN reviewers. Create one listing:
- Test environments (from GitHub Actions matrix)
- R CMD check results (0 errors, 0 warnings, 0 notes)
- Any downstream dependencies

### 14. Create `inst/WORDLIST` for spell checking
sensemakr maintains a custom wordlist (`inst/WORDLIST`) with technical terms and author names that would otherwise be flagged by CRAN spell checks. Create one with terms like: sensemakr, mrsensemakr, Mendelian, pleiotropic, Cinelli, etc.

### 15. Update `.Rbuildignore`
Expand to match sensemakr's pattern. Add:
```
^\.github$
^cran-comments\.md$
^plan\.md$
^CRAN-RELEASE$
^CRAN-SUBMISSION$
^docs$
^_pkgdown\.yml$
^pkgdown$
^LICENSE\.md$
^NEWS\.md$
```

### 16. Update `.gitignore`
sensemakr includes:
```
inst/doc
.DS_Store
.Rproj.user
.Rhistory
.RData
.Ruserdata
tests/testthat/Rplots.pdf
```
Add `tests/testthat/Rplots.pdf` and any other missing entries.

### 17. Regenerate documentation with current roxygen2
`RoxygenNote: 7.1.1` is old. sensemakr uses `7.3.2`. Running `roxygen2::roxygenise()` with a current version will update the NAMESPACE, man pages, and the `RoxygenNote` field.

---

## Optional / Future

### 18. Set up pkgdown documentation site
sensemakr has a full pkgdown site at `http://carloscinelli.com/sensemakr/` hosted via GitHub Pages. Consider setting up the same for mrsensemakr with `_pkgdown.yml` and a GitHub Actions deploy workflow.

### 19. Version bump
Version `0.3` is fine for a first CRAN submission, but consider `0.3.0` following proper semver. sensemakr is at `0.1.6`.

### 20. Run `R CMD check --as-cran`
After all fixes, run the full CRAN check locally to catch any remaining issues before submission.

---

## Implementation Order

**Phase 1 — Package metadata & docs (must fix)**
1. Fix DESCRIPTION (`Authors@R`, Title, Description, URL, BugReports, R version)
2. Add roxygen2 docs for `print.mr_sensemakr` and `plot.mr_sensemakr` (including `@return`)
3. Add `@return` to `mr_sensemakr()` roxygen block
4. Improve `sim_data` documentation (description + variable docs)
5. Fix vignette `T`/`F` usage
6. Fix test file issues

**Phase 2 — CI/CD & infrastructure (matching sensemakr)**
7. Add `.github/workflows/R-CMD-check.yaml`
8. Add `.github/workflows/test-coverage.yaml`
9. Update README badges (remove Travis/AppVeyor, add GHA/Codecov/CRAN)
10. Remove legacy CI files (`.travis.yml`, `appveyor.yml`, `codecov.yml`)

**Phase 3 — Polish**
11. Create `NEWS.md`
12. Create `inst/CITATION`
13. Create `cran-comments.md`
14. Create `inst/WORDLIST`
15. Update `.Rbuildignore` and `.gitignore`
16. Regenerate docs (`roxygen2::roxygenise()`)
17. Run `R CMD check --as-cran` and fix any remaining issues
