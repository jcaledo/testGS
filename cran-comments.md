## Submission

This submission can be considered a first submission, since a previous version of the package (orthGS v0.1.5) was archived. Now, in the current version, the use of the Bioconductor package muscle (previously listed under Suggests) has been discontinued. Additionally, a new vignette has been introduced, providing various alternatives for performing multiple sequence alignments directly within the R environment.

## Previous correspondence related to the submission

Subject: Response to auto-check issues for orthGS version 0.1.6
Dear CRAN Team,
Thank you for your detailed feedback regarding the auto-check results for orthGS version 0.1.6. I have carefully reviewed the noSuggests log and identified the issues related to the examples in the mapTrees function. These issues arose from the use of functionality provided by the fs package, which was listed under Suggests. I have now replaced this functionality with alternative R code that produces similar output without requiring the fs package. A corrected version will be resubmitted promptly.
Regarding the issues reported with the previous CRAN release of orthGS version 0.1.5, they should now be resolved. In this current version, I have addressed the problems by avoiding the dependencies that caused the prior concerns.
Please feel free to let me know if further clarification or adjustments are needed.
Best regards,
Elena Aledo


## Test environments 

* local OS X install, R 4.4.0

* win-builder (devel and release)

* rhubv2: linux (R-devel), macos (R-devel), macos-arm64 (R-devel), windows (R-devel).

## R CMD check results

0 errors | 0 warnings | 0 notes
