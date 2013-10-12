#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
##################################################################
#          Add variable name replacement rules as                #
#                                                                #
#             oldname => newname,                                #
#                                                                #
#          One rule per line!                                    #
#                                                                #
#  This file should be named RenameList.pl and be in the same    #
#  directory as Rename.pl, or it should be explicitly given      #
#  with the -i=FileName switch of Rename.pl                      #
#                                                                #
#                                                                #
##################################################################

%newname=(
#	  ModProcMH   =>   IE_ModProc,
#	  ModMain     =>   IE_ModMain,
#	  ModIO       =>   IE_ModIo,
#	  MPI_COMM_WORLD => iComm,
	  ModIonosphere => IE_ModIonosphere,
	  );
