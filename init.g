#
# The crystallographic groups kit based on GAP related packages and interfaces to other 3rd tools.
#
# Reading the declaration part of the package.
#

ReadPackage( "CrystKit", "gap/CrystKit.gd");

# DeclareAutoreadableVariables failed with "Error, Too many open files (internal file descriptor limit reached)".
# https://mail.google.com/mail/u/0/?ogbl#sent/QgrcJHsbdJTKBzblSRZbvxTpxKrWGxgLfJb
# https://docs.gap-system.org/doc/ref/chap76.html#X8495E5327D563AC3
# ref:
# Public/repo/github.com/gap-system/gap.git/pkg/smallsemi/init.g
DeclareAutoreadableVariables( "CrystKit", "gap/names.g",
                              [ "cryst2names", "cryst3names", "aff3names", "aff4names"  ] );


