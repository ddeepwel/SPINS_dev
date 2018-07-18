# Contributing to SPINS

SPINS is maintained by Chris Subich (@csubich).
If you have made developments to SPINS that you would like to include in the main repository,
first contact Chris to discuss the matter with him.

## Merging Changes

The uWaterloo GitLab system (<https://git.uwaterloo.ca>, the current host of SPINS) supports merge requests.
Once you have settled on changes that you would like to merge, use the following steps to prepare the merge request.
1. Create a fork of the main SPINS repository if you don't already have one
2. Apply your changes to the clean fork. Make sure that you follow good commit standards and
  1. Commit incremental changes when possible
  2. Leave informative commit messages
  3. Update the version number as appropraite (see section on Version Numbers)
3. Once you've commit all of your changes to your fork, go to the main page for SPINS (<https://git.uwaterloo.ca/SPINS/SPINaS_main>)
  1. Select the *Merge Resquests* tab
  2. Choose *New Merge Request*
  3. Select the appropriate fork/branch as the source
  4. Leave an informative merge message that summarizes the changes being merged.


## Version Numbers

Version numbers are tracked in the file VERSION (which are then automatically passed to
the executables at compile-time).
This means that changing the version number is as simple as updating the VERSION file.

When pushing changes to the main SPINS repository, you should increment the version numbers accordingly.

StackOverflow has a nice summary of when to change which version number
(https://stackoverflow.com/questions/3826580/what-rules-does-software-version-numbering-follow)
- Major version numbers change whenever there is some significant change  being introduced. For example, a large or potentially backward-incompatible change to a software package.
- Minor version numbers change when a new, minor feature is introduced or when a set of smaller features is rolled out.
- Patch numbers change when a new build of the software is released to customers. This is normally for small bug-fixes or the like.
