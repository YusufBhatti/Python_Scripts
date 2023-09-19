# script to update multiple Git repositories

 # Array of directories to update
 #declare -a dirs=("Model_Source_Code" "NEMO_trcdms-File" "onedrive" "orca" "Python_Scripts")
 # I want to automatically detect the file names:
 

#for dir in */ ; do
#      echo "Updating $dir ..."
#        cd "$dir" || continue
#          if [ -d .git ]; then
#                  git pull origin main
#                    else
#                            echo "Skipping $dir as it's not a Git repository."
#                              fi
#                                cd ..
#                            done
# Find all directories and loop over them to pull updates
for dir in */ ; do
      echo "Updating $dir ..."
        cd "$dir" || continue
          if [ -d .git ]; then
                  git fetch
                      if [ $? -ne 0 ]; then
                                echo "Failed to fetch updates for $dir. Skipping..."
                                      cd ..
                                            continue
                                                fi
                                                    branch=$(git symbolic-ref --short HEAD)
                                                        git pull origin "$branch"
                                                            if [ $? -ne 0 ]; then
                                                                      echo "Failed to pull updates for $dir on branch $branch."
                                                                          fi
                                                                            else
                                                                                    echo "Skipping $dir as it's not a Git repository."
                                                                                      fi
                                                                                        cd ..
                                                                                    done
