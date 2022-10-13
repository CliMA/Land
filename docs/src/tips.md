# Git and Julia Tips for CliMA.Land


## Rebase the commits before merging into main
1. Switch to the feature branch and make sure you have a clean workspace
   - `git checkout FeatureBranch`
   - `git status` to confirm zero uncommitted changes
2. Backup current branch to avoid unexpected errors
   - `git checkout -b FeatureBranch_Backup`
3. Sync the main branch
   - `git checkout main`
   - `git pull`
4. Switch to the feature branch and merge main into it
   - `git checkout FeatureBranch`
   - `git status` to confirm zero uncommitted changes
   - `git merge main`
   - Resolve conflicts and commit them if any
5. Reabse the feature branch
   - `git reset origin/main` to rewrite history
   - `git diff` (optional) would show all the changes
   - `git add --all` to add local changes on top of main
   - `git commit -m "Single commit message"`
   - `git push -f` to force push because of the re-written history
6. Create a Pull Request
   - Through the web
7. Remove local unnecessary branch (e.g., FeatureBranch_backup)
   - `git branch -d FeatureBranch_Backup`
   - `git branch -D FeatureBranch_Backup` to force remove
