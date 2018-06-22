//##########################################################################
//#                                                                        #
//#                       CLOUDCOMPARE PLUGIN: qPOSC                       #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 of the License.               #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#                             COPYRIGHT: XXX                             #
//#                                                                        #
//##########################################################################



//my
#include "qPOSC.h"
#include "ccPoscDlg.h"

//Qt
#include <QtGui>
#include <QMainWindow>

//qCC_db
#include <ccPointCloud.h>

//system
#include <iostream>

//functions
#include "functions.h"

qPOSC::qPOSC(QObject* parent/*=0*/)
	: QObject(parent)
	, m_action(0)
{
}

void qPOSC::onNewSelection(const ccHObject::Container& selectedEntities)
{
}

void qPOSC::getActions(QActionGroup& group)
{
	if (!m_action)
	{
		m_action = new QAction(getName(),this);
		m_action->setToolTip(getDescription());
		m_action->setIcon(getIcon());
		connect(m_action, SIGNAL(triggered()), this, SLOT(doAction()));
	}

	group.addAction(m_action);
}

void qPOSC::doAction()
{
	assert(m_app);
	if (!m_app)
		return;

    /* window
    ccPoscDlg dlg(m_app->getMainWindow());
    if (!dlg.exec())
        return;
    QString a = QString("[testQQ] %1").arg(dlg.thresholdBox->value());
    m_app->dispToConsole(a,ccMainAppInterface::STD_CONSOLE_MESSAGE);
       window */


    const ccHObject::Container& selectedEntities = m_app->getSelectedEntities();
    size_t selNum = selectedEntities.size();
    if (selNum!=1)
    {
        m_app->dispToConsole("Select only one cloud!",ccMainAppInterface::ERR_CONSOLE_MESSAGE);
        return;
    }
    ccHObject* ent = selectedEntities[0];
    assert(ent);
    if (!ent || !ent->isA(CC_TYPES::POINT_CLOUD))
    {
        m_app->dispToConsole("Select a real point cloud!",ccMainAppInterface::ERR_CONSOLE_MESSAGE);
        return;
    }

    ccPointCloud* pc = static_cast<ccPointCloud*>(ent);

    int pointsNumber = pc->size();
    double **points = newDoubleMatrix(pointsNumber, 3);

    for (int i=0; i<pointsNumber; i++)
    {
        const CCVector3* P = pc->getPoint(i);
        points[i][0] = static_cast<double>(P->x);
        points[i][1] = static_cast<double>(P->y);
        points[i][2] = static_cast<double>(P->z);
    }


    std::cout << "original:\n";
    for(int i = 0; i < pointsNumber; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            std::cout << points[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";


    cov(points, pointsNumber, 3);
    deleteDoubleMatrix(points, pointsNumber);


}

QIcon qPOSC::getIcon() const
{
    return QIcon(":/CC/plugin/qPOSC/icon.png");
}
