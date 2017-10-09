#include "mainwindow.h"
#include "ui_mainwindow.h"



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    for( QCustomPlot* plot: { ui->plotA, ui->plotB})
    {
        plot->addGraph();

        plot->graph(0)->setLineStyle(QCPGraph::lsNone);
        plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 1));

        plot->xAxis->setRange(0, 11);
        plot->yAxis->setRange(0, 11);
    }


    createPoints( 50*0.002 );
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::createPoints(double noise)
{
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<> dist(0.0, noise);

    _points_x.clear();
    _points_y.clear();

    for (int x=1; x<=10; x++)
    {
        for (int y=1; y<=10; y++)
        {
            for (int i=0; i<40; i++)
            {
                _points_x.push_back( x + dist(e2) );
                _points_y.push_back( y + dist(e2) );
            }
        }
    }

    ui->plotA->graph(0)->setData(_points_x, _points_y);
    ui->plotA->replot();

    updateMeanShift();
}

void MainWindow::updateMeanShift()
{
    MeanShift2D::PointsVector points, shifted;
    points.resize(_points_x.size());

    for(int i=0; i<_points_x.size(); i++)
    {
        points[i] = { _points_x[i], _points_y[i] };
    }

    _meanshift.setMeanshiftEps( ui->epsilonSpinbox->value() );
    shifted = _meanshift.meanshift(points, ui->kernelSpinbox->value() );

    _shifted_x.clear();
    _shifted_y.clear();

    for(int i=0; i<shifted.size(); i++)
    {
        _shifted_x.push_back( shifted[i][0] );
        _shifted_y.push_back( shifted[i][1] );
    }
    ui->plotB->graph(0)->setData(_shifted_x, _shifted_y);
    ui->plotB->replot();
}

void MainWindow::on_horizontalSlider_sliderMoved(int position)
{
    double pos = static_cast<double>(position) * 0.002;
    createPoints(pos);
}

void MainWindow::on_epsilonSpinbox_valueChanged(double)
{
    updateMeanShift();
}

void MainWindow::on_kernelSpinbox_valueChanged(double)
{
    updateMeanShift();
}

void MainWindow::on_radioGaussian_toggled(bool checked)
{
    if(checked) {
        _meanshift.setKernelFunction( GaussianKernel );
        updateMeanShift();
    }
}

void MainWindow::on_radioQuartic_toggled(bool checked)
{
    if(checked){
        _meanshift.setKernelFunction( QuarticKernel );
        updateMeanShift();
    }
}

void MainWindow::on_radioParabolic_toggled(bool checked)
{
    if(checked){
        _meanshift.setKernelFunction( ParabolicKernel );
        updateMeanShift();
    }
}
