======================
koopmans documentation
======================

Documentation for the ``koopmans`` code, hosted at https://koopmans-docs.readthedocs.io/en/stable/

Building locally
----------------

In order to build a local copy of the documentation, you will first need to install the required ``python`` packages by running

.. code-block:: bash

   pip install -r requirements.txt

Then to build the documentation simply run

.. code-block:: bash

   make html

If you then open your favourite browser and open ``./_build/index.html`` you will be able to browse the local build of the documentation. N.B. if you subsequently modify and re-build the documentation, you may need to force your browser to reload the page in order to see those changes reflected in the browser.

Contributing
------------

In order to contribute to the website...

1. create a new local branch
2. make your modifications to your local branch
3. push your new branch to ``github``
4. create a pull request on ``github`` into the ``latest`` branch

If your pull request is accepted then your modifications will automatically appear on the `latest <https://koopmans-docs.readthedocs.io/en/latest/>`_ version of the website. N.B. by default the website directs users to the `stable <https://koopmans-docs.readthedocs.io/en/stable/>`_ version of the website; in the bottom left of the website there is an option to switch between versions.

Contact
-------
Written and maintained by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna (2020-)

For help and feedback email edward.linscott@gmail.com

