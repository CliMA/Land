# Git and Julia Tips for CliMA.Land

## Instantiate the project
1. The Land project
   - `cd .../Land` to the home directory of the project
   - `julia --project`
   - `using Pkg`
   - `Pkg.instantiate()`
   - Delete the file `Manifest.toml` when encountering errors
   - Redo `Pkg.instantiate()`

2. The docs/ project
   - `cd .../Land/docs`
   - `julia --project`
   - `using Pkg`
   - `Pkg.instantiate()`
   - Delete the file `Manifest.toml` when encountering erros
   - Delete the line `Land="*"` in the `Project.toml`
   - Type `]` in the julia REPL to enter the pkg environment
   - `dev ..` or `add ..` to add the Land project
   - Redo `Pkg.instantiate()`

## Run local tests before merging into main
1. Optional: load the Land project
   - `cd .../Land` to the home directory of the project
   - `julia --project`
   - `using Land` to make sure the project can be loaded
2. Optional: disable the tutorial to save some time
   - change `generate_tutorials = true` to `generate_tutorials = false` in file `docs/list_of_tutorials.jl`
   - change the line back to `generate_tutorials = true` before commiting the changes
3. Optional: initialize the project
   - `cd .../Land` to the home directory of the project
   - `julia --project -e "using Pkg; Pkg.instantiate(); Pkg.build(); Pkg.precompile();"`
   - `julia --project=docs/ -e "using Pkg; Pkg.instantiate(); Pkg.build(); Pkg.precompile();"`
4. Test the documentation project
   - `julia --project=docs/ docs/make.jl`
   - Fix errors and warnings associated with Land project
5. Test the Land project
   - `julia --project -e "using Pkg; Pkg.test();"`
   - Resolve failed tests

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
