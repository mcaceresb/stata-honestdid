library(fixest)
data(base_stagg)
res_sunab <- feols(y ~ x1 + sunab(year_treated, year) | id + year, base_stagg)

  x <- res_sunab
  sunab_agg   <- x$model_matrix_info$sunab$agg_period
  sunab_names <- names(x$coefficients)
  sunab_sel   <- grepl(sunab_agg, sunab_names, perl=TRUE)
  sunab_names <- sunab_names[sunab_sel]
  if(!is.null(x$weights)){
    sunab_wgt <- colSums(x$weights * sign(model.matrix(x)[, sunab_names, drop=FALSE]))
  } else {
    sunab_wgt <- colSums(sign(model.matrix(x)[, sunab_names, drop=FALSE]))
  }
  sunab_cohorts <- as.numeric(gsub(paste0(".*", sunab_agg, ".*"), "\\2", sunab_names, perl=TRUE))
  sunab_mat     <- model.matrix(~ 0 + factor(sunab_cohorts))
  sunab_trans   <- solve(t(sunab_mat) %*% (sunab_wgt * sunab_mat)) %*% t(sunab_wgt * sunab_mat)
  sunab_coefs   <- sunab_trans %*% cbind(x$coefficients[sunab_sel])
  sunab_vcov    <- sunab_trans %*% x$cov.scaled[sunab_sel, sunab_sel] %*% t(sunab_trans)
  sunab_se      <- sqrt(diag(sunab_vcov))

  xnames <- attr(x$model_matrix_info$sunab$agg, "model_matrix_info")$coef_names_full
  max(abs(x$coeftable[rownames(x$coeftable) %in% xnames,1:2] - cbind(sunab_coefs, sunab_se)))

str(res_sunab)
etable(res_sunab)
sqrt(diag(res_sunab$cov.iid))
sqrt(diag(res_sunab$cov.unscaled))
attr(res_sunab$coeftable, "type")

res_sunab$se
res_sunab$coeftable
res_sunab$sigma2
res_sunab$weights

x <- res_sunab
V = x$cov.scaled
use_weights = TRUE
mm = model.matrix(x)
agg = x$model_matrix_info$sunab$agg_period
is_name = !is.null(names(agg))
coef = x$coefficients
cname = names(coef)
qui = grepl(agg, cname, perl = TRUE)
cname_select = cname[qui]
root = gsub(paste0(".*", agg, ".*"), "\\1", cname_select, perl = TRUE)
val  = gsub(paste0(".*", agg, ".*"), "\\2", cname_select, perl = TRUE)
  name_df = unique(data.frame(root, val, stringsAsFactors = FALSE))
  v_names = cname_select[root == r & val == v]
  c_all = c()
  se_all = c()
  for(i in 1:nrow(name_df)){
    r = name_df[i, 1]
    v = name_df[i, 2]
    v_names = cname_select[root == r & val == v]

    if(use_weights && !is.null(x$weights)){
      shares = colSums(x$weights * sign(mm[, v_names, drop = FALSE]))
    } else {
      shares = colSums(sign(mm[, v_names, drop = FALSE]))
    }

    shares = shares / sum(shares)

    # The coef
    c_value = sum(shares * coef[v_names])

    # The variance
    n = length(v_names)
    s1 = matrix(shares, n, n)
    s2 = matrix(shares, n, n, byrow = TRUE)

    var_value = sum(s1 * s2 * V[v_names, v_names])
    se_value = sqrt(var_value)

    c_all[length(c_all) + 1] = c_value
    se_all[length(se_all) + 1] = se_value
  }




names(summary(res_sunab)$cov.scaled)
names(summary(res_sunab)$se)

names(res_sunab$coefficients)
names(res_sunab$coeftable)
coef(res_sunab)
vcov(res_sunab, "iid")
