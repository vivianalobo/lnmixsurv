# make_survival_mixture <- function() {
#     set_new_model("survival_mixture")
#     set_model_mode(model = "survival_mixture", mode = "censored regression")
#     set_model_engine(
#         "survival_mixture",
#         mode = "censored regression",
#         eng = "ln_mixture"
#     )
#     set_dependency("survival_mixture", eng = "ln_mixture", pkg = "persistencia")

#     survival_mixture <-
#         function(mode = "censored regression") {
#             # Check for correct mode
#             if (mode != "censored regression") {
#                 rlang::abort("`mode` should be 'censored regression'")
#             }

#             # Save some empty slots for future parts of the specification
#             new_model_spec(
#                 "survival_mixture",
#                 args = NULL,
#                 eng_args = NULL,
#                 mode = mode,
#                 method = NULL,
#                 engine = NULL
#             )
#         }

#     set_fit(
#         model = "survival_mixture",
#         eng = "ln_mixture",
#         mode = "censored regression",
#         value = list(
#             interface = "formula",
#             protect = c("formula", "data"),
#             func = c(pkg = "persistencia", fun = "survival_ln_mixture"),
#             defaults = list()
#         )
#     )

#     set_encoding(
#         model = "survival_mixture",
#         eng = "ln_mixture",
#         mode = "censored regression",
#         options = list(
#             predictor_indicators = "traditional",
#             compute_intercept = FALSE,
#             remove_intercept = FALSE,
#             allow_sparse_x = FALSE
#         )
#     )

#     surv_info <-
#         list(
#             pre = NULL,
#             post = NULL,
#             func = c(fun = "predict"),
#             args =
#             # These lists should be of the form:
#             # {predict.mda argument name} = {values provided from parsnip objects}
#                 list(
#                     # We don't want the first two arguments evaluated right now
#                     # since they don't exist yet. `type` is a simple object that
#                     # doesn't need to have its evaluation deferred.
#                     object = quote(object),
#                     new_data = quote(new_data),
#                     type = "survival",
#                     time = quote(time)
#                 )
#         )
#     set_pred(
#         model = "survival_mixture",
#         eng = "ln_mixture",
#         mode = "censored regression",
#         type = "survival",
#         value = surv_info
#     )
# }


# dados <- readRDS("./dados_simulados_ln_mix.rds") %>% mutate(x = factor(ifelse(x == 1, "F", "M")))

# survival_mixture() %>% translate(engine = "ln_mixture")

# sr_spec <-
#     survival_mixture() %>%
#     set_engine("ln_mixture") %>%
#     set_mode("censored regression")
# sr_spec

# set.seed(8975)
# sr_fit <- sr_spec %>% fit(Surv(y, delta) ~ x, dados)
# sr_fit

# options(expressions = 500000)
# predict(sr_fit, new_data = dados[1, ], type = "survival", time = seq(120))
# show_model_info("survival_mixture")
