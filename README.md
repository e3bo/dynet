# dynet
Code from a few research projects pertaining to modeling disease spread on networks

The main dependencies are igraph and the GNU scientific library.

If you would like to use autotools to build it, run the autogen.sh
script. Then you can

    ./configure
    make
    make check
    make install

as usual. Please be forewarned that in spite of the polished
appearance of some parts of this program, it is experimental code and
likely crawling with bugs.

Also note that the uthash.h file is a third-party component that is
separately copyrighted and licensed.


