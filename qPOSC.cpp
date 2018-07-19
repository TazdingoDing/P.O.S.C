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
#include "SPG.h"

//Qt
#include <QtGui>
#include <QMainWindow>

//qCC_db
#include <ccPointCloud.h>
#include <ccPolyline.h>

//system
#include <iostream>

//functions
#include "functions.h"


/*	QCC_DB_LIB_API extern const Rgba white;
    QCC_DB_LIB_API extern const Rgba lightGrey;
    QCC_DB_LIB_API extern const Rgba darkGrey;
    QCC_DB_LIB_API extern const Rgba red;
    QCC_DB_LIB_API extern const Rgba green;
    QCC_DB_LIB_API extern const Rgba blue;
    QCC_DB_LIB_API extern const Rgba darkBlue;
    QCC_DB_LIB_API extern const Rgba magenta;
    QCC_DB_LIB_API extern const Rgba cyan;
    QCC_DB_LIB_API extern const Rgba orange;
    QCC_DB_LIB_API extern const Rgba black;
    QCC_DB_LIB_API extern const Rgba yellow;
 */
static ccColor::Rgb colors[3] = {ccColor::red,ccColor::black,ccColor::blue };
static int current = 0;

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
    if (selNum<1)
    {
        m_app->dispToConsole("Select at least one cloud!",ccMainAppInterface::ERR_CONSOLE_MESSAGE);
        return;
    }



    for (int i = 0; i < selNum; i++)
    {
        ccHObject* ent = selectedEntities[i];
        assert(ent);
        if (!ent || !ent->isA(CC_TYPES::POINT_CLOUD))
        {
            //m_app->dispToConsole("Select a real point cloud!",ccMainAppInterface::ERR_CONSOLE_MESSAGE);
            return;
        }



        ccPointCloud* pc = static_cast<ccPointCloud*>(ent);
        int pointsNumber = pc->size();
        if (pointsNumber%2 != 0)
        {
            m_app->dispToConsole("points number not 2x",ccMainAppInterface::ERR_CONSOLE_MESSAGE);
            return;
        }

        ccColor::Rgb color = colors[current++];
        current %= 3;

        ccPolyline **poly = new ccPolyline*[pointsNumber/2];
        for (int j =0; j < pointsNumber; j+=2)
        {
            int index = j/2;
            poly[index] = new ccPolyline(pc);
            poly[index]-> addPointIndex(j,j+2);
            poly[index]->setColor(color) ;
            poly[index]->setWidth(4);
            pc->addChild(poly[index]);
            m_app->addToDB(poly[index]);
            poly[index]->showColors(true);

        }
        pc->redrawDisplay();
    }





    /*

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

    ccPointCloud *newCloud = getNewCloud(points);
*/




/*   color
    if (pc->reserveTheRGBTable())
    {
        pc->setRGBColor(static_cast<ColorCompType>(255),static_cast<ColorCompType>(0),static_cast<ColorCompType>(0) );
        m_app->setSelectedInDB(ent, false);
        pc->showSF(false);
        pc->showColors(true);
        pc->redrawDisplay();
        //m_app->redrawAll();
        //m_app->setSelectedInDB(ent, true);
    }
color   */


/*  PCA
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
PCA  */




    //m_app->doActionAddConstantSF();




    //deleteDoubleMatrix(points, pointsNumber);


}

QIcon qPOSC::getIcon() const
{
    return QIcon(":/CC/plugin/qPOSC/icon.png");
}
