

use gnuplot::*;

pub fn figure(low_res_time: &[f64], low_res_brownian: &[f64],
              high_res_time: &[f64], high_res_brownian: &[f64]) {

    let col1 = "#000000";
    let col2 = "#2CA02C";

    let mut fg = Figure::new();
    let ax = fg.axes2d();

    ax.set_x_grid(true);
    ax.set_y_grid(true);

    ax.lines(low_res_time, low_res_brownian, &[Color(col1)]);
    ax.points(low_res_time, low_res_brownian, &[Color(col1), PointSymbol('o')]);


    ax.lines(high_res_time, high_res_brownian, &[Color(col2)]);
    ax.points(high_res_time, high_res_brownian, &[Color(col2), PointSymbol('o')]);



    // set smoother gnuplot terminal and show it
    fg.set_terminal("wxt", "");

    // spawn plot on new thread to let Rust code continue
    std::thread::spawn(move || {
        fg.show().unwrap();
    });

}
