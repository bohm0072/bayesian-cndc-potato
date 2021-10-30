# Init ------------------------

library(usethis)
library(gitcreds)

# Create GitHub token ---------

# Use GitHub authetication and create personal access token [PAT] via browser interface
usethis::create_github_token(
  scopes = c("repo", "user", "gist", "workflow"),
  description = "R:GITHUB_PAT",
  host = NULL
)

# Using R terminal interface, set PAT created from authentication above
gitcreds::gitcreds_set()

# END -------------------------