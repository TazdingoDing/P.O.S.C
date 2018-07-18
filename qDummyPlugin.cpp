//##########################################################################
//#                                                                        #
//#                       CLOUDCOMPARE PLUGIN: qDummy                      #
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

//First: replace all occurrences of 'qDummyPlugin' by your own plugin class name in this file!
#include "qDummyPlugin.h"

//Qt
#include <QtGui>

//qCC_db
#include <ccPointCloud.h>
#include <ccPolyline.h>
#include <ccBBox.h>
#include <cc2DLabel.h>

//system
#include <iostream>
#include <string>
#include <limits>

//Default constructor: should mainly be used to initialize
//actions (pointers) and other members
qDummyPlugin::qDummyPlugin(QObject* parent/*=0*/)
	: QObject(parent)
	, m_action(0)
{
}

//This method should enable or disable each plugin action
//depending on the currently selected entities ('selectedEntities').
//For example: if none of the selected entities is a cloud, and your
//plugin deals only with clouds, call 'm_action->setEnabled(false)'
void qDummyPlugin::onNewSelection(const ccHObject::Container& selectedEntities)
{
	//if (m_action)
	//	m_action->setEnabled(!selectedEntities.empty());
}

//This method returns all 'actions' of your plugin.
//It will be called only once, when plugin is loaded.
void qDummyPlugin::getActions(QActionGroup& group)
{
	//default action (if it has not been already created, it's the moment to do it)
	if (!m_action)
	{
		//here we use the default plugin name, description and icon,
		//but each action can have its own!
		m_action = new QAction(getName(),this);
		m_action->setToolTip(getDescription());
		m_action->setIcon(getIcon());
		//connect appropriate signal
		connect(m_action, SIGNAL(triggered()), this, SLOT(doAction()));
	}

	group.addAction(m_action);
}

float* getMinMaxXYZ(ccPointCloud* pc)
{
    unsigned pointsNum = pc->size();
    float x, y, z, minX, maxX, minY, maxY, minZ, maxZ;
    minX = minY = minZ = std::numeric_limits<float>::max();
    maxX = maxY = maxZ = std::numeric_limits<float>::min();

    for (unsigned i=0; i<pointsNum; ++i)
    {
        const CCVector3* P = pc->getPoint(i);
        x = static_cast<float>(P->x);
        y = static_cast<float>(P->y);
        z = static_cast<float>(P->z);
        if (x > maxX)
        {
            maxX = x;
        }
        if (x < minX)
        {
            minX = x;
        }
        if (y > maxY)
        {
            maxY = y;
        }
        if (y < minY)
        {
            minY = y;
        }
        if (z > maxZ)
        {
            maxZ = z;
        }
        if (z < minZ)
        {
            minZ = z;
        }
    }

    float* minMaxXYZ = new float[6];
    minMaxXYZ[0] = minX;
    minMaxXYZ[1] = maxX;
    minMaxXYZ[2] = minY;
    minMaxXYZ[3] = maxY;
    minMaxXYZ[4] = minZ;
    minMaxXYZ[5] = maxZ;

    return minMaxXYZ;
    //std::cout<< "minx = "<< minX << ", miny = "<< minY <<", minz = "<< minZ <<"\n";
    //std::cout<< "maxx = "<< maxX << ", maxy = "<< maxY <<", maxz = "<< maxZ <<"\n";
}


ccPolyline* getPolylineBox(double minX, double maxX,double minY, double maxY,double minZ, double maxZ)
{
    ccPointCloud *boxPointCloud = new ccPointCloud("8 Points");

    CCVector3 *points = new CCVector3[8];
    points[0] = CCVector3(minX, minY, maxZ);
    points[1] = CCVector3(minX, maxY, maxZ);
    points[2] = CCVector3(maxX, maxY, maxZ);
    points[3] = CCVector3(maxX, minY, maxZ);
    points[4] = CCVector3(minX, minY, minZ);
    points[5] = CCVector3(minX, maxY, minZ);
    points[6] = CCVector3(maxX, maxY, minZ);
    points[7] = CCVector3(maxX, minY, minZ);

    boxPointCloud->reserve(8) ;
    for (int i=0; i<8; i++)
    {
        boxPointCloud->addPoint(points[i]);
    }

    ccPolyline *poly = new ccPolyline(boxPointCloud) ;

    poly->addPointIndex(0,4);
    poly->addPointIndex(0);
    poly->addPointIndex(4,6);
    poly->addPointIndex(1);
    poly->addPointIndex(5,7);
    poly->addPointIndex(2);
    poly->addPointIndex(6,8);
    poly->addPointIndex(3);
    poly->addPointIndex(7);
    poly->addPointIndex(4);

    ccColor::Rgb red = ccColor::red ;
    poly->setColor( red ) ;

    return poly;
}



//This is an example of an action's slot called when the corresponding action
//is triggered (i.e. the corresponding icon or menu entry is clicked in CC's
//main interface). You can access to most of CC components (database,
//3D views, console, etc.) via the 'm_app' attribute (ccMainAppInterface
//object).
void qDummyPlugin::doAction()
{
	//m_app should have already been initialized by CC when plugin is loaded!
	//(--> pure internal check)
	assert(m_app);
	if (!m_app)
		return;

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

    ccBBox boxx = pc->getOwnBB();

    CCVector3 minCorner = boxx.minCorner();
    CCVector3 maxCorner = boxx.maxCorner();
    float minX, maxX, minY, maxY, minZ, maxZ;
    minX = static_cast<float>(minCorner.x);
    minY = static_cast<float>(minCorner.y);
    minZ = static_cast<float>(minCorner.z);
    maxX = static_cast<float>(maxCorner.x);
    maxY = static_cast<float>(maxCorner.y);
    maxZ = static_cast<float>(maxCorner.z);
    std::cout<< "minx = "<< minX << ", miny = "<< minY <<", minz = "<< minZ <<"\n";
    std::cout<< "maxx = "<< maxX << ", maxy = "<< maxY <<", maxz = "<< maxZ <<"\n";

/*
    ccBBox bbox = pc->getBoundingBox();
    CCVector3 minCorner = bbox->minCorner();
    CCVector3 maxCorner = bbox->maxCorner();

    float minX, maxX, minY, maxY, minZ, maxZ;
    minX = static_cast<float>(minCorner->x);
    minY = static_cast<float>(minCorner->y);
    minZ = static_cast<float>(minCorner->z);
    maxX = static_cast<float>(maxCorner->x);
    maxY = static_cast<float>(maxCorner->y);
    maxZ = static_cast<float>(maxCorner->z);
    std::cout<< "minx = "<< minX << ", miny = "<< minY <<", minz = "<< minZ <<"\n";
    std::cout<< "maxx = "<< maxX << ", maxy = "<< maxY <<", maxz = "<< maxZ <<"\n";
*/
    float* minMaxXYZ = getMinMaxXYZ(pc);

    ccPolyline *box = getPolylineBox(minMaxXYZ[0],minMaxXYZ[1],minMaxXYZ[2],minMaxXYZ[3],minMaxXYZ[4],minMaxXYZ[5]);

    pc->addChild(box);
    m_app->addToDB(box);

    std::cout<< "done" <<"\n";

    QString title = QString("QQ");
    cc2DLabel *label = new cc2DLabel();
    box->addChild(label);
    m_app->addToDB(label);
    label->setVisible(pc->isVisible());
    label->setDisplayedIn2D(false);
    label->addPoint(pc, pc->size() - 1);
    label->setName(title);
    label->displayPointLegend(true);


/*  print points
    unsigned count = pc->size();
    float x, y, z;
    for (unsigned i=0; i<count; ++i)
    {
        m_app->dispToConsole("[qDummyPlugin] point",ccMainAppInterface::STD_CONSOLE_MESSAGE);
        const CCVector3* P = pc->getPoint(i);
        x = static_cast<float>(P->x);
        y = static_cast<float>(P->y);
        z = static_cast<float>(P->z);
        std::string s_x = std::to_string(x);
        std::string s_y = std::to_string(y);
        std::string s_z = std::to_string(z);
        std::cout<< "x = "<< s_x << ", y = "<< s_y <<", z = "<<s_z <<"\n";
        //m_app->dispToConsole(a,ccMainAppInterface::STD_CONSOLE_MESSAGE);
    }

    std::cout<< count <<"\n";
*/


/*  draw box
    double minX, maxX, minY,maxY, minZ, maxZ;
    minX=200;
    maxX=230;
    minY=160;
    maxY=220;
    minZ=20;
    maxZ=30;

    ccPolyline *box = getPolylineBox(minX,maxX,minY,maxY,minZ,maxZ);

    pc->addChild(box);
    m_app->addToDB(box);
*/


    //m_app->dispToConsole("[qDummyPlugin] HeQQ world!",ccMainAppInterface::STD_CONSOLE_MESSAGE); //a standard message is displayed in the console
    //m_app->dispToConsole("[qDummyPlugin] Warning: dummy plugin shouldn't be used as is!",ccMainAppInterface::WRN_CONSOLE_MESSAGE); //a warning message is displayed in the console
    //m_app->dispToConsole("Dummy plugin shouldn't be used as is!",ccMainAppInterface::ERR_CONSOLE_MESSAGE); //an error message is displayed in the console AND an error box will pop-up!

	/*** HERE ENDS THE ACTION ***/

}





//This method should return the plugin icon (it will be used mainly
//if your plugin as several actions in which case CC will create a
//dedicated sub-menu entry with this icon.
QIcon qDummyPlugin::getIcon() const
{
	//open qDummyPlugin.qrc (text file), update the "prefix" and the
	//icon(s) filename(s). Then save it with the right name (yourPlugin.qrc).
	//(eventually, remove the original qDummyPlugin.qrc file!)
	return QIcon(":/CC/plugin/qDummyPlugin/icon.png");
}
