# bioGUItools
Tools for bioGUI (e.g. CWL Parser/argpase2biogui/...)

bioGIU Tools is a collection of simple helpers for anyone wanting to create new templates for bioGUI.
The content of this repository is described below.

## x2biogui

An python argparse-extension which first supplies two useful argparse-type extensions and second provides means to convert the python argparse parser into a bioGUI template.

### FileStubType
Similar to argparse.FileType(...) but does not specify the full file, but the location to a (or possibly many) file(s) as well as the file's prefix.

E.g. when storing files into the directory results, and each file should have the same prefix, but multiple analyses were performed.

The program would take the FileStub-Location, e.g. results/condition1_ and creates, e.g. the files results/condition1_{detail.png|summary.png}
                      
### FolderType

Similar to the argparse.FileType(...), but accepts as input any valid folder.
This is essential when the user should specify a path.

### cwl2biogui

This tool is based on the cwl2argparse tool provided by [Common Workflow Language](https://github.com/common-workflow-language/cwl2argparse).
