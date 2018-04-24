require "rubygems"
require "tmpdir"

require "bundler/setup"
require "jekyll"


# Change this as needed
GITHUB_REPONAME = "dnanexus/dnanexus-example-applets"


namespace :site do
  desc "Generate GitHub page files"
  task :generate do
    Jekyll::Site.new(Jekyll.configuration({
      "source"      => ".",
      "destination" => "_site",
      "baseurl"     => ""
    })).process
  end


  desc "Generate and publish blog to gh-pages"
  task :publish => [:generate] do
    puts "Verifying git tree is clean"
    verify_git_tree_clean()
    Dir.mktmpdir do |tmp|
      cp_r "_site/.", tmp

      pwd = Dir.pwd
      Dir.chdir tmp

      system "touch .nojekyll"  # Prevent Github from trying to render based on _folders
      system "git init"
      system "git add -A"  # note -A used, untracked files are added
      message = "Site updated at #{Time.now.utc}"
      system "git commit -m #{message.inspect}"
      system "git remote add origin git@github.com:#{GITHUB_REPONAME}.git"
      system "git push origin master:refs/heads/gh-pages --force"

      Dir.chdir pwd
    end
  end
end

def verify_git_tree_clean()
  # TODO verify on master and up to date with origin master
  verified = ''
  # Check working
  verified += `set -e; git diff-files --quiet || (echo 'Working directory has changes'; exit 1)`

  # Check staging
  verified += `set -e; git diff-index --quiet --cached HEAD -- || (echo 'There are staged, but not commited files'; exit 1)`

  # Untracked and unignored files
  verified += `set -e; test -z "$(git ls-files .. --exclude-standard --others)" || (echo "There are untracked files"; exit 1)`

  if verified != ''
    puts verified
    fail
  end
end