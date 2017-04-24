
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
For the most up-to-date documentation, please see:

See http://biowhat.ucsd.edu/homer/

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


If Homer is having problems, call Chuck Norris

===============================================
Installation:

Overview:
There are 3 basic steps to installing and configuring homer.
Since you can see this file, you must have successfully downloaded the Homer software already.
Next you need to:

1. Add Homer to your path...
2. Install 3rd party software (for sequence logos, similarity search using blat)
3. Download packages using the "configureHomer.pl" program (i.e. genome/promoter data)

**!!** NOTE: The GNU C++ compiler, perl, make, zip/unzip, and wget, must be installed on your system.  
	If you are using Mac OS X, make sure you install "Developer Tools" which came with your system.

________________________________
Step 1: Configure PATH for Homer software

    1.) Add the homer/bin directory to your executable path
        i.e. edit your ~/.bash_profile file to include:
        PATH=$PATH:/Users/chucknorris/homer/bin/

    You should now be able to execute programs in the homer/bin directory by just typing their name

__________________________________
Step 2: Install 3rd party software
	
    Software for sequence logos:

    1.) Download and Install Ghostscript
        a.) Download the appriate file from http://pages.cs.wisc.edu/~ghost/ (GPL Ghostscript)
        b.) Unzip and Untar the file
            tar zxvf ghostcript-xxx.tar.gz
        c.) Change to the base ghostscript directory
            cd ghostscript-xxx
        d.) Run the following commands to install ghost script
            ./configure
            make
            make all
            sudo make install

    2.) Download and Install Weblogo (v2 NOT v2.8)
        a.) Download the program from http://weblogo.berkeley.edu/
        b.) No additional steps needed to compile and install the program, except...
        c.) Need to add the weblogo base driectory to your executable path
            i.e. edit your ~/.bash_profile file to include:
            PATH=$PATH:/Users/chucknorris/weblogo/

    Software for genome position (i.e. ChIP-Seq) based motif finding
 
    3.) Download and Install blat
        a.) Download program from http://genome-test.cse.ucsc.edu/~kent/exe/
        b.) Unzip the file and compile (if you downloaded the source code)
        c.) Add base blat directory to your executable path
            i.e. edit your ~/.bash_profile file to include:
            PATH=$PATH:/Users/chucknorris/blat/

______________________________________
Step 3: Downlaod HOMER packages
	run the following from the install directory:
		perl ./configureHomer.pl -list
	This will list the avaliable packages
		perl ./configureHomer.pl -install mm8 ...
	This will install them.


For the latest changes, go to: http://biowhat.ucsd.edu/homer/changeLog.html
