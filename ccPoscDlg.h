//##########################################################################
//#                                                                        #
//#                       CLOUDCOMPARE PLUGIN: qPOSC                       #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 or later of the License.      #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#                  COPYRIGHT: xxx                                        #
//#                                                                        #
//##########################################################################

#ifndef CC_POSC_DLG_HEADER
#define CC_POSC_DLG_HEADER

#include <ui_poscDlg.h>

//! Dialog for the PCV plugin
class ccPoscDlg : public QDialog, public Ui::POSCDialog
{
public:
	explicit ccPoscDlg(QWidget* parent = 0);
};

#endif
