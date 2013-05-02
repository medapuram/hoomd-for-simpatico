# Maintainer: joaander

################################
## Version information
set(HOOMD_VERSION_MAJOR "0")
set(HOOMD_VERSION_MINOR "11")
set(HOOMD_VERSION_PATCH "3")
set(HOOMD_VERSION "${HOOMD_VERSION_MAJOR}.${HOOMD_VERSION_MINOR}.${HOOMD_VERSION_PATCH}")

# users may not have git installed, or this may be a tarball build - set a dummy version if that is the case
include(GetGitRevisionDescription)
git_describe(HOOMD_GIT_VERSION)
if (HOOMD_GIT_VERSION)
    set(HOOMD_VERSION_LONG "${HOOMD_GIT_VERSION}")
else (HOOMD_GIT_VERSION)
    set(HOOMD_VERSION_LONG "${HOOMD_VERSION}-unknown")
endif (HOOMD_GIT_VERSION)

get_git_head_revision(GIT_REFSPEC GIT_SHA1)
if (GIT_REFSPEC)
    set(HOOMD_GIT_REFSPEC "${GIT_REFSPEC}")
else (GIT_REFSPEC)
    set(HOOMD_GIT_REFSPEC "unknown")
endif (GIT_REFSPEC)

if (GIT_SHA1)
    set(HOOMD_GIT_SHA1 "${GIT_SHA1}")
else (GIT_SHA1)
    set(HOOMD_GIT_SHA1 "unknown")
endif (GIT_SHA1)

