#include <QApplication>
#include "interface.h"

int main(int argc,char** argv)
{
    QTextCodec::setCodecForLocale(QTextCodec::codecForName("UTF-8"));
    
    QApplication app(argc,argv);
    MainWindow *o = new MainWindow;
    return app.exec();
}
