# Setting up GitHub access from JASMIN

After getting access to the hk26 GitHub repository by following instructions at [python_setup_conda](https://github.com/digital-earths-UK-hackathon/hk26/blob/main/docs/python_setup_conda.md) you may wish to push changes back to the hk26 or other relevant GitHub repository.

To do this you will need to be setup with GitHub access from Jasim.

1. Generate a new ssh key and add it to your GitHub credentials, by following instructions at
[GitHub pages](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent?platform=linux).

```
ssh-keygen -t ed25519 -C "your_email@example.com"
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

Then add to your GitHub profile online, following instructions [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).

2. Update your `~/.ssh/config` file with following text:
```
Host github.com
    Hostname ssh.github.com
    Port 443
    User git
```

3. If setting up GitHub access for the first time, you may need to also type in terminal the following 1-time commands:
```
git config --global user.name "Your GitHub User Name"
git config --global user.email "your-email@example.com"
```

4. Test you can check connection to GitHub with:
```
ssh -T git@github.com
```


To enable pushing local changes to the hk26 repository, you will also need to

5. Check the path of the remote repo:
```
git remote -v
```
If SSH is configured, this will return
```
origin	git@github.com:digital-earths-UK-hackathon/hk26.git (fetch)
origin	git@github.com:digital-earths-UK-hackathon/hk26.git (push)
```

Else, if HTTPS is configured it may give
```
origin  https://github.com/digital-earths-UK-hackathon/hk26.git (fetch)
origin  https://github.com/digital-earths-UK-hackathon/hk26.git (push)
```
In this case, you will need to type:
```
git remote set-url origin git@github.com:digital-earths-UK-hackathon/hk26.git
```

6. See e.g. [CSET documentation](https://metoffice.github.io/CSET/contributing/git.html#useful-git-commands) for a summary of useful git commands.

In general if making local changes to a file, you will need to action the following 3 steps to ADD the modified file to next commit, COMMIT changes to local repository, and then PUSH the change to the hk26 repo to share with others:
```
git add <filename>
git commit
git push
```
