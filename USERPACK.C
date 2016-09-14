#include "d:\work\bc\include\mcadincl.h"
    
    extern FUNCTIONINFO FDNA;

    
    char *ErrorMessageTable[] = {
    "argument must be real",    //  error 1 --  argument must be real
    "insufficient memory",      //  error 2 --  memory allocation error
    "interrupted",              //  error 3 --  execution interrupted
    "Ошибка открытия файла ДНА" //  error 4 -- ошибка в функции FReadDNA_KVAFunction
    };


BOOL WINAPI DllEntryPoint (HANDLE hDLL, DWORD dwReason, LPVOID lpReserved)
{
  switch (dwReason)
  {
    case DLL_PROCESS_ATTACH:
    {

      // DLL is attaching to the address space of the current process.
      //
		  if (!CreateUserErrorMessageTable( hDLL, 4, ErrorMessageTable ) )
				break;

 		  if ( CreateUserFunction( hDLL, &FDNA) == NULL )
 				break;
        

    }

	 case DLL_THREAD_ATTACH:        // A new thread is being created in the current process.
	 case DLL_THREAD_DETACH:        // A thread is exiting cleanly.
	 case DLL_PROCESS_DETACH:      // The calling process is detaching the DLL from its address space.

		break;
  }
  return TRUE;
}



