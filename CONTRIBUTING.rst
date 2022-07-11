Guidelines for contributing to koopmans
=======================================

The ``koopmans`` project is `hosted on GitHub <https://github.com/epfl-theos/koopmans>`_ as a public repository, owned the *epfl-theos* GitHub organization.

To contribute to the project:

1. create a branch
2. make changes to your branch
3. push these changes to GitHub
4. create a pull request

Your proposed changes will be reviewed by the developers and, if accepted, merged into the master branch

Signing up for GitHub
---------------------

To contribute, you first need a GitHub account

* If you do not already have a GitHub account, `sign up here <https://github.com/join>`_
* If you already have a GitHub account, we recommend that you get an `education account <https://education.github.com/benefits?type=student>`_ if you are eligible
  
In order to connect to GitHub using SSH (recommended), you will need to `add an SSH key to your account <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_

To contribute to the koopmans repository, you need to become a collaborator. Please contact one of the team administrators with your GitHub username to request this. Current administrators include:

* Edward Linscott
* Riccardo De Gennaro
* Nicola Colonna

Once you are given collaborator status you will have write access to the official repository, allowing you to create branches. Collaborators cannot directly push changes to the master branch, but must contribute these via pull requests.

Creating a pull request
-----------------------

When you are working on a contribution, open a pull request early on in the process. It does not need to be ready to be reviewed: simply `mark it as a draft <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/changing-the-stage-of-a-pull-request#converting-a-pull-request-to-a-draft>`_. Having a draft pull request open is useful because it will give you access to useful tools such as automated testing and coding style suggestions. It is also helpful for us as developers to keep track of who is working on what.

Once your pull request is ready, mark your request as `ready for review <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/changing-the-stage-of-a-pull-request#marking-a-pull-request-as-ready-for-review>`_ and request a review from one of the developers. If they accept your pull request, it will be merged into the master branch.

If there are issues with your pull request, you may be asked to make changes. There is no need to create a new pull request. You can apply the requested changes directly to the branch for which the pull request was issued and they will be automatically added your pull request.

Using pre-commit
^^^^^^^^^^^^^^^^
We recommend you use the tool ``pre-commit`` in order to perform some automated checks before you make commits. To set this up, simply run ``pre-commit install``.

Working on the Quantum ESPRESSO submodule
-----------------------------------------
By default, the official Quantum ESPRESSO repository is a submodule of ``koopmans``, located in the subdirectory ``quantum_espresso/q-e``. Because we do not administer the Quantum ESPRESSO repository, hotfixes to bugs will first be published to `this fork <https://gitlab.com/elinscott/qe_koopmans>`_ of Quantum ESPRESSO.

If you want to obtain a hotfix, add this fork as a remote to your local repository by running ``git remote add qe_koopmans git@gitlab.com:elinscott/qe_koopmans.git`` from within the ``q-e`` subdirectory. You will then be able to pull changes from this fork to your local repository via ``git pull qe_koopmans``.

N.B. If you have write-access to this fork, you can also push changes to it via ``git push qe_koopmans <branch>``.

