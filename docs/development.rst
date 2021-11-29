.. _Development:

===================================
Contributing to Further Development
===================================

Nemo is hosted on `github <https://github.com/simonsobs/nemo/>`_ and is available under a free
software license. Help is appreciated in its development.


Reporting Issues
----------------

Found a bug? Please tell us, either on the `issues page <https://github.com/simonsobs/nemo/issues>`_
in the Simons Observatory `#sz Slack channel <https://simonsobs.slack.com/messages/C35CDSEGJ>`_
(if you are an SO member), or by `email <hiltonm@ukzn.ac.za>`_.


Contributing Code
-----------------

Want to add a feature or fix something? There are some tasks that need attention listed on the 
`issues page <https://github.com/simonsobs/nemo/issues>`_.

The preferred method of contributing code is to clone the `repository <https://github.com/simonsobs/nemo>`_ 
and work on your new feature in your own branch::

    git clone https://github.com/simonsobs/nemo.git
    git checkout -b the-name-of-your-branch

When the time comes to commit your changes, please contact `Matt Hilton <hiltonm@ukzn.ac.za>`_  with your
github username, in order to be granted write access to the repository. You only need to do this once.

The ``main`` branch itself is protected. When you are ready for your changes to be added to it, please
issue a Pull Request for your changes to be reviewed. 


Style
^^^^^

When adding code, please adhere to the style used throughout Nemo where possible (see below).

* As you may notice, Nemo uses `camelCase <https://en.wikipedia.org/wiki/Camel_case>`_ throughout -
  please keep it that way.

* Indent with 4 spaces.

* The maximum line length is 110 characters (sometimes it makes sense to break this).

* Docstrings *should* follow the `Google style <https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html>`_.
  This has only been partially done (so far) in the existing code. At the very least, every function
  should have a docstring of some kind that describes what it does, even if it is re-formatted later.


Testing
^^^^^^^

Nemo uses the Robot framework for tests (see :ref:`TestingPage`). These are integration tests,
rather than unit tests (at least at the moment), and can take a while to run. While more work 
(and more tests) need to be added, you should check that these tests still pass (or at least, 
don't crash) before committing your changes.


Checking Memory Usage
^^^^^^^^^^^^^^^^^^^^^

If you need to check the memory usage of Nemo, we recommend the use of the 
`memory-profiler <https://pypi.org/project/memory-profiler/>`_ package. For example,
to use it on Nemo::
    
    mprof run /usr/local/bin/nemo test_MFMF_S16_auto.yml

This will log memory usage over time, which can be plotted using other tools in the 
`memory-profiler <https://pypi.org/project/memory-profiler/>`_ package.
