- Zero Y values are not removed for clustering
- For variable selection, weights are used in cv.glmnet
- `type.measure` is set equal to `mae` in cv.glmnet
- After selecting features, rows of X are not subsetted for fitting regression