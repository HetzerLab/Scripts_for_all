# Quick steps for using git and github

# 1. Create a GitHub account, then send me or swati you github ID so we can add you to the lab's GitHub page.

# 2. Create new repository on GitHub. To avoid errors, do not initialize the new repository with README, license, or gitignore files. 

# 3. In the terminal, change the current working directory to your local project. Initialize the local directory as a Git repository.

git init

# 4. Add the files in your new local repository. This stages them for the first commit.

git add .

# 5. Commit the files that you've staged in your local repository.

git commit -m "First commit"

# 6. At the top of your GitHub repository's Quick Setup page, click to copy the remote repository URL. In Terminal, add the URL for the remote repository where your local repository will be pushed.

git remote add origin remote repository URL

git remote -v
# Verifies the new remote URL

# 7. Push the changes in your local repository to GitHub.

git push -u origin master
