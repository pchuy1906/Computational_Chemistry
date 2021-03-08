git init 
git add *
git commit -m "Hand on RDKit"

git branch -M master

git remote set-url origin https://github.com/pchuy1906/Computational_Chemistry.git
git remote add     origin https://github.com/pchuy1906/Computational_Chemistry.git

git fetch origin master
git merge origin master

git push --force origin master


##…or create a new repository on the command line
##echo "# Data_Science_Learning" >> README.md
##git init
##git add README.md
##git commit -m "first commit"
##git branch -M main
##git remote add origin https://github.com/pchuy1906/Data_Science_Learning.git
##git push -u origin main
##
##…or push an existing repository from the command line
##git remote add origin https://github.com/pchuy1906/Data_Science_Learning.git
##git branch -M main
##git push -u origin main
