# Software & Dependency Guide

This document provides direct links and permission instructions for the external tools required by the Automated Docking Pipeline.

## 1. Official Software Links
* **AutoDock Vina**: [Official GitHub Releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
    * *Download the version specific to your OS (e.g., `vina_1.2.5_mac_x86_64` or `vina_1.2.5_linux_x86_64`).*
* **UCSF ChimeraX**: [Download Page](https://www.cgl.ucsf.edu/chimerax/download.html)
* **OpenBabel**: [Official Documentation](http://openbabel.org/wiki/Main_Page)

---

## 2. Creating the ChimeraX Symlink
The script looks for a local executable named `chimerax`. You must create a symbolic link (symlink) in the project folder pointing to your installed application.

### **macOS**
Run this in your terminal (update the version number to match your installation):
```bash
ln -s /Applications/ChimeraX-1.8.app/Contents/MacOS/ChimeraX ./chimerax
```

### **linux**
```bash
ln -s /usr/bin/chimerax ./chimerax
```

## 3. AutoDock Vina: Setup and Permissions

After downloading the AutoDock Vina executable, follow these steps to prepare it for the pipeline.

### **Step A: Rename the Executable**
The script specifically looks for a file named `vina`. Rename your downloaded binary based on your OS:

**macOS/Linux**: 
  ```bash
  mv vina_1.2.5_mac_arm64 vina
  ```
### **Step B: Grant Execution Permission**

 **macOS/Linux**: 
  ```bash
 chmod +x vina
  ```

  ### **Step C: Verification**

 **macOS/Linux**: 
  ```bash
 ./vina --version
  ```
