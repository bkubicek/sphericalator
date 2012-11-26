#include "window.h"
#include "glwidget.h"
#include "sphericalcalculator.h"
#include "patch.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <iostream>
#include <QPushButton>
#include <QString>
#include <QtGui>

Window::Window()
{
    QHBoxLayout *hb1=new QHBoxLayout();
    //sb_n=new QSpinBox();

    saving=false;
    glw=new GLWidget(this);
    QVBoxLayout *vb=new QVBoxLayout();
    ovs.resize(5);
    for(int i=0;i<ovs.size();i++)
    {
        ovs[i]=new OverlayGUI();
        if(i==0)
        {
            ovs[i]->intensity->setValue(1);
            ovs[i]->l->setDisabled(true);
            ovs[i]->m->setDisabled(true);
        }
        connect(ovs[i],SIGNAL(changed()),this,SLOT(updateOverlays()));
        vb->addLayout(ovs[i]);
    }

    pb_export=new QPushButton("export");
    vb->addWidget(pb_export);
    vb->addStretch();
    connect(pb_export,SIGNAL(clicked()),this,SLOT(exportFile()));
    //vb->addWidget(sb_n);

    hb1->addLayout(vb);
    hb1->addWidget(glw);

    QWidget *central=new QWidget();
    central->setLayout(hb1);
    setCentralWidget(central);
    sc= new SphericalCalculator();
    updateOverlays();
    recalculate();
    //connect(sb_n,SIGNAL(valueChanged(int)),this,SLOT(recalculate(int)));

}
void Window::updateOverlays()
{
    sc->overlays.resize(ovs.size());
    std::cout<<"Recalc Overlays"<<std::endl;
    for(int i=0;i<ovs.size();i++)
    {
        ovs[i]->toOverlay(sc->overlays[i]);
    }
    recalculate();
}

void Window::recalculate(int i)
{
    recalculate();
}

void Window::recalculate()
{

    int maxl=0;
    for(int i=0;i<ovs.size();i++)
        if(ovs[i]->l->value()>maxl)
            maxl=ovs[i]->l->value();
    if(maxl<4)
        maxl=4;
    int maxm=0;
    for(int i=0;i<ovs.size();i++)
        if(ovs[i]->m->value()>maxl)
            maxm=ovs[i]->m->value();
    if(maxm<4)
        maxm=4;

    if(!saving)
        sc->setResolution(maxl*20,maxm*35);
    else
        sc->setResolution(maxl*35,maxm*45);
    saving=false;

    sc->update();

    glw->glp->clear();
    //std::cout<<"Faces:"<<sc->faces.size()<<std::endl;
    for(int i=0;i<(int)sc->faces.size();i++)
    {
        Face &f=sc->faces[i];
        QVector3D a(f.x[0][0],f.x[0][1],f.x[0][2]);
        QVector3D b(f.x[1][0],f.x[1][1],f.x[1][2]);
        QVector3D c(f.x[2][0],f.x[2][1],f.x[2][2]);

        QVector3D n;
        n=QVector3D::crossProduct(b-a,c-a);
        float r=0.2;
        glw->glp->addTri(r*a,r*b,r*c,-r*n);
    }

    //glw->glp->addTri(QVector3D(0,0,0),QVector3D(2,2,2),QVector3D(5,5,2),QVector3D(-1,0,0));
    glw->glp->setSmoothing(Patch::Smooth);

    glw->g->finalize();
    glw->updateGL();
}

void Window::exportFile()
{
    QString filename=QFileDialog::getSaveFileName(this, "Save file", "*.stl", "*.stl");
    saving=true;
    recalculate();
    sc->saveFile(filename.toStdString().c_str());

}

OverlayGUI::OverlayGUI():QHBoxLayout()
{

    l=new QSpinBox();
    l->setMaximum(19);
    m=new QSpinBox();
    m->setMinimum(0);
    intensity=new QDoubleSpinBox();
    intensity->setRange(-1,5);
    intensity->setSingleStep(0.05);
    addWidget(new QLabel("l:"));
    addWidget(l);
    addWidget(new QLabel("m:"));
    addWidget(m);
    addWidget(new QLabel("Amplitude:"));
    addWidget(intensity);
    connect(l,SIGNAL(valueChanged(int)),this,SLOT(change(int)));
    connect(m,SIGNAL(valueChanged(int)),this,SLOT(change(int)));
    connect(intensity,SIGNAL(valueChanged(double)),this,SLOT(change(double)));
    setStretch(0,0);
}
void OverlayGUI::change()
{
    m->setMaximum(l->value());

    std::cout<<"Change emitted"<<std::endl;
    emit changed();
}

void OverlayGUI::change(int i)
{
    change();
}
void OverlayGUI::change(double i)
{
    change();
}

void OverlayGUI::toOverlay(Overlay &o)
{
    o.intensity=intensity->value();
    o.l=l->value();
    o.m=m->value();
}

