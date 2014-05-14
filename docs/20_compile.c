/**
   \page page20  Compile the Code

   This step can be skipped if the user choose to use the precompiled binary
   files, which are available for all the popular plateforms: Windows, Mac, and
   Linux. Please notice that the job scheduling and monitoring utilities are not
   included in the released binaries.

   \section sect-requirement Requirements

    - C99 compliant compiler: GCC4, Intel ICC, or llvm-clang.

    - (Optional) GTK+ and Cairo. For drawdaemon and monitor.

    - (Optional in Linux/Mac) FFTW version 3. Will download from MAOS site if not available in the system.
    
    - (Optional in Linux/Mac) Optimized blas and lapack. Blas and lapack can
    usually be installed through the linux distribution official
    repository. Will download from MAOS site if not available in the system.
    
    
    \section sect-prep Preparing the folders and compiling

    We recommend using three different folders to 1) store the source tree, 2)
    compile the code, and 3) run the simulation.

    The source tree is where the source is extracted. Obtain the software in the
    form of the compressed tar ball: maos_version.tar.gz. Extract it into some
    directory and name the directory src_dir. Example:

    \verbatim
    cd ~/work/programming
    tar zxvf maos_version.tar.gz 
    cd maos_version
    export src_dir=$(pwd)
    \endverbatim
    
    Next we create another folder, where we are going to compile the code. example:

    \verbatim
    cd ~/work/maos
    mkdir comp_optim && cd comp_optim
    $src_dir/configure
    make -j4 #this will compile using 4 threads.
    \endverbatim

    which will config the compiling folder with GCC compiler and default
    optimization flags and compile the code.
    
    To use Intel C++ compiler,

    \verbatim
    cd ~/work/maos
    $src_dir/configure CC=icc
    make
    \endverbatim
    
    To for using the redistributed Intel MKL library (for fast cholesky decomposition)
    \verbatim
    $src_dir/configure --enable-mkl
    \endverbatim

    To enable debug mode
    \verbatim
    $src_dir/configure --enable-debug
    \endverbatim
    
    To specify your own compiler using gcc4.1 which is not the system default
    \verbatim
    CC=gcc4.1 $src_dir/configure
    \endverbatim

    For clusters that use job management tools, disable the built-in scheduler by using
    \verbatim
    $src_dir/configure --disable-scheduler
    \endverbatim
    
    The compiled executable is maos in the sub-folder “bin” of the compiling
    folder. You do not have to do "make install" to run the simulations.

    \subsection sect-mac-gtk Installing GTK+ in MAC OS and Compile Monitor, Drawdaemon

    First, if xcode is not installed in full, make sure the following packages are installed
    \verbatim
    DeveloperToolsCLI.pkg
    gcc4.2.pkg
    MACOSX10.5.pkg
    DevSDK.pkg
    OpenGLSDK.pkg
    DeveloperToolsSystemSupport
    \endverbatim
    Second, download and install pkg-config from pkgconfig.freedesktop.org/releases/     

    Third, following the instructions on page
    gtk-osx.sourceforge.net/ to install gtk+ for osx. The 
    environmental variables are set by jhbuild, here we take important
    ones and put them in our .bash_profile to avoid running jhbuild
    each time. Put the following in your .bash_profile. Adjust the path accordingly if you changed gtk-osx installing options.
    \verbatim
    export PKG_CONFIG_PATH="/usr/lib/pkgconfig:$HOME/gtk/inst/lib/pkgconfig:$HOME/gtk/inst/share/pkgconfig"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/gtk/inst/lib"
    export PATH="$PATH:$HOME/gtk/inst/bin"
    export C_INCLUDE_PATH="$HOME/gtk/inst/include"
    \endverbatim
    Now rerun MAOS configure and make. Monitor and drawdaemon should appear in bin folder now.

    To make the app bundles. 
    First, install quartz engine and ige-mac-buildler. 
    sourceforge.net/apps/trac/gtk-osx/wiki/Bundle 
    jhbuild build gtk-quartz-engine
    
    Second, symbolic link bin/drawdaemon and bin/monitor to $HOME/gtk/inst/bin
    And cd to scripts/bundle folder, run build.sh. .app files should appear in your desktop.
    
    \subsection sect-monitor Using the job monitor
    
    When reasonably recent GTK+ and Cairo libraries are present in the system,
    two additional executives will be compiled, “drawdaemon” and “monitor”. The
    plotting utility “drawdaemon” is located in sub-folder “lib”. It is launched
    automatically when plotting commands are issued in the software. The
    monitoring utility “monitor” is located in sub-folder “bin”.  Link monitor to
    a folder that is in your PATH and type “monitor” in the command line to use
    it.

    In order to use monitor to monitor jobs in different machines, follow these steps:

    Create a folder .aos in your home folder, create two files in there. Put
    these files in all the machines you want to monitor.
    \verbatim
~/.aos/hosts  Each line contains a machine's hostname you want to monitor. 
              They should be the names returned by command "hostname", 
              and you should be able to log into those machines 
              with these names. Set up /etc/hosts properly to ensure this.
~/.aos/port   contains the port you want the scheduler to bind to 
              so that the monitor can connect to it. 
              Any non-previledged port numbers are fine (like 10000). 
              Different users should have different port number. 
	      The port has to be the same on all the machines. 
    \endverbatim
*/
