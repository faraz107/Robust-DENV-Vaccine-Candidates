
export_formattable <- function(f, file, zoom, width, height, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay, 
          zoom = zoom)
}

custom_color_tile <- function (...) 
{
  formatter("span",
            style = function(x) style(display = "block", 
                                      padding = "0 4px", 
                                      `color` = "black", 
                                      font.size = "1.5em",
                                      `border-radius` = "4px", 
                                      `background-color` = csscolor(gradient(as.numeric(x), 
                                                                             ...))))
}

custom_title <- function (...) 
{
  formatter("span",
            style = function(x) style(`color` = "black", 
                                      font.size = "1.5em"))
}
