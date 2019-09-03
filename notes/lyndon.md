Lyndon White:ox:  16 minutes ago
Least squares linear regression + suitable transformation then.
f(Y) = wX solve for w
Where Y and X are matrixes containing observations of your variates/co-variates.
The vector w is a weighted feature mapping.

Kevin Bonham:microbiome:  10 minutes ago
Cool - I'll give this a shot. Thanks!

Lyndon White:ox:  3 minutes ago
There may well be formulations of this which constraint w to be nonnegative (if there are not, then you can code one up in JuMP, or just work with abs(w) and fit it with Flux or optim).
And assuming negative contributions are not allowed (which might require transforming your inputs and outputs).
that would be much more interpretable.
And it would be related to nonnegative matrix factorisation. Which is related to LDA
