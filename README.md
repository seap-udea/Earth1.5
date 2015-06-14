# Earth1.5
Global Environmental Model including Biotic Regulation

Useful links
------------

Presentation
------------

This is Earth1.5.

Getting a copy
--------------

To get a copy of the newest version of this project just execute:

```
$ git clone --branch master http://github.com/seap-udea/plynet.git
```

If you want to get a different branch of the project just change
"master" by the name of the branch.

This is the developer copy.  You may also obtain the latest release of
the installable package, available in the `dist` folder.

Instructions for the contributor
--------------------------------

1. Generate a public key of your account at the client where you will
   develop contributions:
   
   ```
   $ ssh-keygen -t rsa -C "user@email"
   ```

2. Upload public key to the github Seap-Udea repository (only authorized
   for the Seap-Udea repository manager), https://github.com/seap-udea.

3. Configure git at the client:

   ```
   $ git config --global user.name "Your Name"
   $ git config --global user.email "your@email"
   ```

4. Get an authorized clone of the project:

   ```
   $ git clone git@github.com:seap-udea/plynet.git
   ```

5. [Optional] Checkout the branch you are interested in
   (e.g. <branch>):

   ```
   $ git checkout -b <branch> origin/<branch>
   ```

6. Checkout back into the master:

   ```
   $ git checkout master
   ```

Licensing
---------

This project must be used and distributed under the [GPL License
Version 2] (http://www.gnu.org/licenses/gpl-2.0.html).

The symbol `[)]` means that it has been developed under the principles
of the [copyleft philosophy](http://en.wikipedia.org/wiki/Copyleft).

All wrongs reserved to [Jorge
I. Zuluaga](mailto:jorge.zuluaga@udea.edu.co).
