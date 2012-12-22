#ifndef WINDOW_H
#define WINDOW_H
#include <QWidget>
#include <QMainWindow>
#include <QHBoxLayout>
class GLWidget;

#include "sphericalcalculator.h"
class Patch;

class QHboxLayout;
class QVboxLayout;
class QSpinBox;
class QCheckBox;
class QDoubleSpinBox;
class QPushButton;

class OverlayGUI:public QHBoxLayout
{
    Q_OBJECT
public:
    OverlayGUI();
    QSpinBox *l;
    QSpinBox *m;
    QDoubleSpinBox *intensity;

    void toOverlay(Overlay &o);
protected:
    inline QSize maximumSize() const{return QSize(200,30);};
signals:
        void changed();
public slots:
    void change();
    void change(int);
    void change(double);



};

class Window:public QMainWindow
{
    Q_OBJECT

public:
    Window();
private:
   GLWidget *glw;

   SphericalCalculator *sc;

   QSpinBox *sb_n;
   std::vector<OverlayGUI*> ovs;
   QPushButton *pb_export;
   
   QCheckBox *support;
   QCheckBox *flat;
   QDoubleSpinBox *supportgap;
   QDoubleSpinBox *flatheight;
   QDoubleSpinBox *diameter;
   
   QCheckBox *autoResolution;
   QSpinBox *meshResolv;
   
   bool saving;
public slots:
   void recalculate();
   void recalculate(int i);
   void recalculate(double d);
   void recalculate(bool b);
    void updateOverlays();
    void exportFile();
};

#endif // WINDOW_H
