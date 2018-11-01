# How to login to the server:

ssh yoursalkid@tardis.salk.edu
## Add password when requested and on the first login agree to trust the certificate.
### Bonus points, you can create an .ssh/config file in your computer to simplify this login process if you'd like.
___________________________________________________________________________________________________________________

# You should save all your data and analysis results in the following folder:

/data/users_data/Yournamehere
## The other folders are mainly for installing things.

### IMPORTANT: The server is not backed up, so make sure you use it for your analysis but save your files somewhere else.

## You can find example scripts to help you get started in:

/data/users_data/Scripts_for_all
____________________________________________________________________________________________________________________

# How to edit you .bash_profile so you can add new folders to your path:

vi .bash_profile

## Then click I to edit the text, when done, click Esc to stop editing.
## Use :save .bash_profile to save your changes, then :q! to leave the text editor

## We have already added a few things to your path, now let's try to add your own data folder here and save.
## You have:

PATH=$PATH:$HOME/.local/bin:$HOME/bin:/bin/STAR/source/:/data/python36/anaconda3:/data/python36/anaconda3/bin:/data:/data/python36:/data/genomes/STARgenomes:/opt/bin:/opt/bedtools2:/opt/bedtools2/bin:/opt/FastQC:/opt/scripts:/opt/tools:/opt/tools/ucsctools:/opt/tools/vcftools:/opt/tools/sratoolkit.2.9.2-centos_linux64:/opt/tools/sratoolkit.2.9.2-centos_linux64/bin

## To add a new folder to this path, so the server knows to go looking for things you want to use in there, add : followed by the path to the folder of interest, like:

:/data/users_data/Juliana

## So that now you have:

PATH=$PATH:$HOME/.local/bin:$HOME/bin:/bin/STAR/source/:/data/python36/anaconda3:/data/python36/anaconda3/bin:/data:/data/python36:/data/genomes/STARgenomes:/opt/bin:/opt/bedtools2:/opt/bedtools2/bin:/opt/FastQC:/opt/scripts:/opt/tools:/opt/tools/ucsctools:/opt/tools/vcftools:/opt/tools/sratoolkit.2.9.2-centos_linux64:/opt/tools/sratoolkit.2.9.2-centos_linux64/bin:/data/users_data/Juliana

### Bonus points, you can use something other than vi if you'd like to, for example nano or vscode. But you need to install and configure those so they will work with the server.

__________________________________________________________________________________________________________________________________
# Tips:

## You can navigate to folders of interest using: 
cd foldername

## To find out the complete path for a folder you are currently in (in case you want to add that to path for example) use 
pwd

## To create new folders use
mkdir foldername

## To delete files use 
rm filename

## If it's a folder use 
rm -rf foldername

## To move files
mv -v /your/from/folder/* /your/to/folder
### The * above means all files in the folder, you can replace that with a specific file name or a different regex

## To search for things:
find / -name "thething"

## To copy a file already in one server folder to another
cp filename /your/to/folder/

## To copy files from a website
wget yoururlhere.com






