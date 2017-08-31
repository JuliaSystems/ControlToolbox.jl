# Contributing to ControlToolbox

First, we would like to thank you for taking the time to contribute to
`ControlToolbox`.

Below we will list, in ascending order of difficulty, some of the ways you can
help us make `ControlToolbox` a better package.

## Bug Reports and Feature Requests

Probably the easiest way to contribute to `ControlToolbox` is through reporting a
bug or a lack of a feature you have noticed when using the package.

Before reporting a bug or requesting a new feature, please note that we try to
collect them in our [Kanban][cc-kanban]. If you have checked the Kanban and could
not find any report related to your problem, please go ahead and [submit][new-issue]
a new issue.

We have provided an issue template for convenience. You do not need to stick to
the template, but we would really appreciate if you filed the issue by giving a
short **description** of the problem, how you think the **expected behaviour**
should look like (hopefully with a small `Julia` code excerpt), how the **actual
behaviour** looks like in the current implementation as well as some **proposals**
and/or **references**, if you can find any, to help solve the issue.

[cc-kanban]: https://github.com/JuliaSystems/ControlToolbox.jl/projects/1
[new-issue]: https://github.com/JuliaSystems/ControlToolbox.jl/issues/new

## Pull Requests

If you would like to get your hands dirty by enhancing our documentation, proposing
new tutorials or walkthroughs, or fixing bugs or implementing new features, you
are more than welcome! We would be more than happy to see you on board :)

Currently, we are following the [conversational development][conv-devel] practice,
which is rather common in `git`-based repositories:

1.  We [submit][new-issue] a new issue for our work, if it is not already
    [listed][cc-kanban]. This way, we inform each other about our contribution,
2.  We [fork][cc-fork] the repository (in)to our own username space,
3.  We create a new branch having a name `xxx-yyy-zzz` where `xxx` is the issue
    number assigned and `yyy-zzz` is some short description,
4.  We do the changes in our own fork and test the changes locally in `Julia`,
5.  We [create][cc-pull] a pull request and elaborate on our solution.

This helps us not only focus our attention on issues which are not being dealt
with by another contributor, but also learn from each other when discussing under
the corresponding pull request.

Please note that we also provide a pull request template which is similar to the
one below:
```
### Changes

Short summary of the changes:

- Change 1. Intuition behind it,
- Change 2. Intuition behind it,
- ...

### Reference

Fixes #issue_number.
```

Since the commits in `git` already show the differences that are made, it would
help if we summarized the changes in simple words together with the intuitions
behind them. Using the `Fixes #issue_number.` in the end is useful to close the
corresponding issue automatically upon a successful merge.

[conv-devel]: https://youtu.be/iV7mVGPXrxU?t=16m25s
[cc-fork]: https://github.com/JuliaSystems/ControlToolbox.jl/fork
[cc-pull]: https://github.com/JuliaSystems/ControlToolbox.jl/pull/new/master

### Documentation Enhancements

A relatively easy way to contribute via a pull request is through writing missing
documentation elements or improving already available ones.

Mainly, we have two different places for our documentation in `ControlToolbox`:

1.  Code documentation in `Julia`. This is what `Julia` uses to show help to users
    issuing a `?command_name`. We try to follow the documentation
    [guidelines][julia-doc] written by the `Julia` team, and we write them in the
    corresponding `.jl` file in the repository,
2.  Documentation in [readthedocs][cc-rtfd]. We are using [reStructuredText][rst-doc]
    for this purpose and keep all the documentation related files under the `doc/`
    folder.

If you choose to contribute to `ControlToolbox` through documentation enhancements,
you can follow steps 1, 2, 3 and 5 above, and you can write a `[ci-skip]`
statement in your commit message to disable testing.

[julia-doc]: http://docs.julialang.org/en/latest/manual/documentation/
[cc-rtfd]: http://controltoolbox.rtfd.io/
[rst-doc]: http://docutils.sourceforge.net/docs/user/rst/quickref.html

### Proposing Tutorials and Walkthroughs

Instead of enhancing documentation in `ControlToolbox`, you might be willing to
create a tutorial showcasing the use of `Julia` and `ControlToolbox` in control
systems identification and analysis applications. In this case, you are more than
welcome to contribute to our [tutorials repository][js-ctj].

[js-ctj]: https://github.com/JuliaSystems/CTJ.git
