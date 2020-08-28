Arnoldi Iteration

This is a demonstration of what the "Arnoldi iteration" looks like when
programmed in ForTran. It is meant to help those who look at the Wikipedia
page, for example, and get lost in the indices. NO fear kids; be grateful
that it is just math and there is nothing ambiguous about it. I wrote it
in ForTran because the parentheses for two-dimensional arrays make it much
easier to understand. For those who see this for the first time, it is a
way to make "nice" vectors to help work with an operator/matrix, possibly
for solving a system of equations, etc. The "nice" hides two good things:
(a) the fact that we use the operator/matrix itself to form the vectors,
and therefore the vectors have an intimate relationship with the operator,
and (b) the Modified Gram-Schmidt makes those vectors form an orthonormal
basis, which is good for sooooo many reasons. If you find this useful, send
chocolate.

Here is the Wikipedia page: https://en.wikipedia.org/wiki/Arnoldi_iteration

Yes, I will delete this repository when I am seek of looking at it. No, I
will not write it in C for you.

IN 2020/08/27

