
#include <QtGui>
#include <QApplication>
#include "mainwindow.h"

using namespace std;

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    MainWindow w;
    w.show();

    int result = app.exec();

    return result;
}
