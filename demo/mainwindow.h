#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <array>
#include "MeanShift/MeanShift.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_horizontalSlider_sliderMoved(int position);

    void updateMeanShift();

    void on_epsilonSpinbox_valueChanged(double arg1);

    void on_kernelSpinbox_valueChanged(double arg1);

    void on_radioGaussian_toggled(bool checked);

    void on_radioQuartic_toggled(bool checked);

    void on_radioParabolic_toggled(bool checked);

private:
    Ui::MainWindow *ui;

    QVector<double> _points_x,  _points_y;
    QVector<double> _shifted_x, _shifted_y;

    void createPoints(double noise);

    typedef MeanShift<std::array<double,2>> MeanShift2D;

    MeanShift2D _meanshift;

};


template <> inline
double euclidean_distance_sqr(const std::array<double,2>& point_a, const std::array<double,2>& point_b){

    const double d0 = (point_a[0] - point_b[0]);
    const double d1 = (point_a[1] - point_b[1]);
    return d0*d0 + d1*d1;
}

template <> inline
void reset_point(std::array<double,2>& point)
{
    point[0] = point[1] = 0.0;
}

template <> inline
void multiply_and_accumulate(const std::array<double,2>& point, double multiplier,  std::array<double,2>& result)
{
    result[0] += point[0] * multiplier;
    result[1] += point[1] * multiplier;
}

template <> inline
void multiply(const std::array<double,2>& point, double multiplier,  std::array<double,2>& result)
{
    result[0] = point[0] * multiplier;
    result[1] = point[1] * multiplier;
}

template <> inline
void add(const std::array<double,2>& point_A, const std::array<double,2>& point_B,  std::array<double,2>& result)
{
    result[0] = point_A[0] + point_B[0];
    result[1] = point_A[1] + point_B[1];
}


#endif // MAINWINDOW_H
