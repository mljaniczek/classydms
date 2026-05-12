# Declare global variables that R CMD check would otherwise flag:
#
# `..` is torch's tensor-indexing placeholder used in expressions like
#   `x[i, ..]` (equivalent to `x[i, , , ]` for a 4-D tensor).
# `self` is the receiver inside `torch::nn_module(initialize=, forward=)`
#   closures, analogous to `this` in other languages. Both are bound by
#   the torch runtime, not by ordinary R scoping, so R CMD check cannot
#   see them statically.
utils::globalVariables(c("..", "self"))
