# fixedtimeanalyzR

This package implements several methods for comparing survival curves at a fixed point in time. The tests are Chi-Square-Tests that rely on test statistics developed by John P. Klein et al. in their article "Analyzing survival curves at a fixed point in time" published in Statistics in Medicine 26 (2007).

## Installation instruction

Currently, there are three ways in order to install **fixedtimeanalyzR**.

### Installation using devtools via GitHub

**Step 1:** Make sure that you have access to the IMBI Heidelberg repository `github.com/imbi-heidelberg`.

**Step 2:** If you haven't already installed `devtools`, install it via
```{r, eval=FALSE}
install.packages("devtools")
```

**Step 3:** In the GitHub settings, you need to generate an access token via **Settings > Developer settings > Personal access tokens > Generate new token**. Don't forget to tick `repo` in order to get access to the private repository that fixedtimeanalyzR is stored in. Copy your access token to clipboard.

**Step 4:**  Install **fixedtimeanalyzR** via
```{r, eval=FALSE}
devtools::install_github(repo="https://github.com/imbi-heidelberg/fixedtimeanalyzR",
                         auth_token="<your personal access token>",
                         build_vignettes=TRUE)
```
We assume that you also want to build the package vignettes. If this is not the case, you may omit the option `build_vignettes=TRUE`.

Success! You may now load the package via `library(fixetimeanalyzR)`.

### Installation using devtools and the package directory

**Step 1:** Download or copy the whole package directory `fixedtimeanalyzR\` to your hard drive.

**Step 2:** If you haven't already installed `devtools`, install it via
```{r, eval=FALSE}
install.packages("devtools")
```

**Step 3:** In a shell, switch to the package directory. Run an R console and type
```{r, eval=FALSE}
devtools::install(build_vignettes=TRUE)
```
We assume that you also want to build the package vignettes. If this is not the case, you may omit the option `build_vignettes=TRUE`.

Success! You may now load the package via `library(fixetimeanalyzR)`.

### Installation using RStudio and the .tar.gz file

**Step 1:** Download or copy the .tar.gz file to your hard drive.

**Step 2:** Launch RStudio and switch to the **Packages** tab. Click **Install > Install from: Package Archive File (.zip, .tar.gz)** and select the .tar.gz file.

Success! You may now load the package via `library(fixetimeanalyzR)`.

## Browsing the vignettes

You may the browse the vignettes by running:
```{r}
browseVignettes("fixedtimeanalyzR")
```
A small tutorial can be found in the introduction vignette.
