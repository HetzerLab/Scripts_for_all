# Quick steps for using git and github

# 1. Create a GitHub account, then send me or swati you github ID so we can add you to the lab's GitHub page.

# 2. Create new repository on GitHub. To avoid errors, do not initialize the new repository with README, license, or gitignore files. 

# 3. In the terminal, change the current working directory to your local project. Initialize the local directory as a Git repository.

git init

# 4. Add the files in your new local repository. This stages them for the first commit.

git add .

# 5. Commit the files that you've staged in your local repository.

git commit -m "New commit 2"

# 6. At the top of your GitHub repository's Quick Setup page, click to copy the remote repository URL. In Terminal, add the URL for the remote repository where your local repository will be pushed.

git remote add origin https://github.com/HetzerLab/Scripts_for_all.git

# 7. Push the changes in your local repository to GitHub.

git push -u origin master

# Repeat steps 5 and 7 whenever you make changes to your files and want to send them to github. 
# If you create new files you want to send over, repeat from 4, 5 and 7.

# If you or someone else made changes to a file uploaded that to github, you should grab the latest changes. For that use:

git fetch origin master

### Cool github cheatsheet: https://education.github.com/git-cheat-sheet-education.pdf

# If you'd like to copy a repository that is on github to your computer:

git clone https://github.com/HetzerLab/Scripts_for_all.git


