remark.macros.upper = function () {
  // `this` is the value in the parenthesis, or undefined if left out
  return this.toUpperCase();
};

remark.macros.scale = function (percentage) {
  var url = this;
  return '<img src="' + url + '" style="width: ' + percentage + '" />';
};

remark.macros.color = function (color) {
  return '<span style="color:' + color + '">' + this + '</span>'
}

remark.macros.emoji = function (emoji) {
  return '<i class="em em-' + emoji + '"></i>'
}

